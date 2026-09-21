"""Bound the validation set by chain length so the validation pass fits in HBM.

Training crops every chain to the preset's crop size (256 residues), but
OpenFold's validation path does not crop: it runs the full chain. Template
triangular attention is roughly cubic in length, so a single long chain decides
the peak memory of the entire validation pass. On this split one 4,174-residue
chain asked for 98.31 GiB on an 80 GB H100 and took the run down with it --
after six hours of training, at the first validation, having written no
checkpoint.

The fix is to exclude the tail rather than to crop it. Cropping would change
what the reported lDDT means; excluding keeps the metric a plain full-chain
lDDT over a slightly smaller, explicitly recorded set.

Only the chain-data cache is rewritten. OpenFoldSingleDataset enumerates
validation chains from the alignment directory and then drops any chain absent
from the cache (data_modules.py, the `c in self.chain_data_cache` filter), so
pruning the cache is sufficient and leaves every symlink in place. The original
cache is kept alongside, so the cap can be raised or lifted without rebuilding
anything.

Read-only with respect to the release, the mmCIF tree and the MSA store.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)

    return digest.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", required=True, type=Path,
                    help="val_mmcif_data_cache.json to prune in place")
    ap.add_argument("--max-length", type=int, default=768,
                    help="longest validation chain to keep, in residues")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    original = json.loads(args.cache.read_text())
    lengths = {k: len(v["seq"]) for k, v in original.items()}

    keep = {k: v for k, v in original.items()
            if lengths[k] <= args.max_length}
    drop = sorted((lengths[k], k) for k in original if k not in keep)

    print(f"cache          {args.cache}")
    print(f"cap            {args.max_length} residues")
    print(f"before         {len(original):,} chains "
          f"(longest {max(lengths.values()):,})")
    print(f"after          {len(keep):,} chains "
          f"(longest {max(lengths[k] for k in keep):,})")
    print(f"excluded       {len(drop):,} chains "
          f"({100 * len(drop) / len(original):.2f}%)")

    if drop:
        shown = ", ".join(f"{k}({n})" for n, k in drop[-10:])
        print(f"longest excluded: {shown}")

    if not keep:
        raise SystemExit("cap excludes every validation chain; refusing")

    if args.dry_run:
        print("\ndry run -- nothing written")
        return 0

    backup = args.cache.with_suffix(".json.uncapped")

    # Keep exactly one pristine copy. Re-running with a different cap must not
    # overwrite the original with an already-pruned one.
    if not backup.exists():
        shutil.copy2(args.cache, backup)
        print(f"\noriginal kept  {backup.name}")

    args.cache.write_text(json.dumps(keep))

    manifest = {
        "stage": "validation_length_cap",
        "max_length_residues": args.max_length,
        "chains_before": len(original),
        "chains_after": len(keep),
        "chains_excluded": len(drop),
        "longest_kept": max(lengths[k] for k in keep),
        "longest_excluded": drop[-1][0] if drop else None,
        "excluded_chains": [k for _, k in drop],
        "uncapped_cache": backup.name,
        "uncapped_cache_sha256": sha256_file(backup),
        "capped_cache_sha256": sha256_file(args.cache),
        "reason": (
            "OpenFold validation does not crop; template triangular attention "
            "is ~cubic in chain length and the tail exhausted an 80 GB H100."
        ),
    }
    manifest_path = args.cache.with_name("validation_length_cap.json")
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True))

    print(f"capped cache   {args.cache.name}")
    print(f"manifest       {manifest_path.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
