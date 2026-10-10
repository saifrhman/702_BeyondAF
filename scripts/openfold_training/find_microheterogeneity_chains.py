"""Find training chains whose projection cannot reproduce the sequence its MSA used.

OpenFold builds a chain's seqres by walking `_entity_poly_seq`. An entry with
point microheterogeneity carries more than one row for the same residue number
-- 9ixd entity 1 models residue 86 three ways (CYS / CSO / CSD) -- and OpenFold
appends every row. Each extra row shifts all later indices, so projecting onto
the Gold-retained residue numbers reads the wrong residues from that point on.
The length can still look right, which is why the input-view adapter compares a
digest rather than a count.

Training then dies the moment the sampler happens to draw such a chain, which
on this run meant 12-16 hours of H100 time per occurrence. Finding them all up
front costs one pass over the mmCIF headers.

Two stages, cheap then exact:

  1. Scan `_entity_poly_seq` of every entry a training chain uses and flag the
     entities carrying duplicate residue numbers. This is a text scan of one
     loop -- no structure is built.
  2. For flagged chains only, run the real projection and compare the digest
     the adapter would compare, so the verdict is the adapter's own, not a
     proxy for it.

Writes the confirmed-bad chain ids so they can be excluded from the training
filter. Read-only with respect to the release, the mmCIF tree and the MSA store.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path


def duplicated_entities(cif_path: Path) -> set[str]:
    """Entity ids in this entry whose _entity_poly_seq repeats a residue number."""
    import gemmi

    try:
        block = gemmi.cif.read(str(cif_path)).sole_block()
    except Exception:
        return set()

    seen: dict[tuple[str, str], int] = {}
    dup: set[str] = set()

    for row in block.find("_entity_poly_seq.", ["entity_id", "num"]):
        key = (row[0], row[1])
        seen[key] = seen.get(key, 0) + 1
        if seen[key] > 1:
            dup.add(row[0])

    return dup


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--projection-index", required=True, type=Path)
    ap.add_argument("--mmcif-dir", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path,
                    help="JSON report: flagged and confirmed-bad chains")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--num-shards", type=int, default=1)
    ap.add_argument("--confirm", action="store_true",
                    help="run the real projection on flagged chains")
    # The two stages need different interpreters: stage 1 needs gemmi, which
    # lives in the PDBClean env, and stage 2 needs Biopython via
    # openfold.data.mmcif_parsing, which lives in the OpenFold env. Neither env
    # has both, so the scan writes its flagged list and the confirm pass reads
    # it back rather than the two sharing a process.
    ap.add_argument("--flagged-input", type=Path, default=None,
                    help="skip the scan and confirm the flagged chains in this "
                         "report (written by a previous --output)")
    args = ap.parse_args()

    index = json.loads(args.projection_index.read_text())

    # group chains by entry so each mmCIF is read once
    by_entry: dict[str, list[str]] = {}
    for chain in index:
        by_entry.setdefault(chain.split("_")[0], []).append(chain)

    entries = sorted(by_entry)
    flagged: list[str] = []
    unreadable: list[str] = []
    mine: list[str] = []

    if args.flagged_input is not None:
        prior = json.loads(args.flagged_input.read_text())
        flagged = list(prior["flagged_chains"])
        unreadable = list(prior.get("unreadable_entries", []))
        print(f"flagged chains loaded from {args.flagged_input.name}: "
              f"{len(flagged):,}", flush=True)
    else:
        mine = [e for i, e in enumerate(entries)
                if i % args.num_shards == args.shard]

        print(f"entries total   : {len(entries):,}")
        print(f"this shard      : {len(mine):,} "
              f"({args.shard + 1} of {args.num_shards})", flush=True)

        for i, pdb in enumerate(mine, 1):
            path = args.mmcif_dir / f"{pdb}.cif"
            if not path.is_file():
                unreadable.append(pdb)
                continue

            dup = duplicated_entities(path)
            if dup:
                for chain in by_entry[pdb]:
                    if index[chain].get("e") in dup:
                        flagged.append(chain)

            if i % 2000 == 0:
                print(f"  {i:,} / {len(mine):,} entries   flagged so far: "
                      f"{len(flagged)}", flush=True)

        print(f"\nflagged (entity has duplicate residue numbers): {len(flagged):,}")

    confirmed: list[dict] = []

    if args.confirm and flagged:
        sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
        from pdbclean.openfold_training_view import (  # noqa: E402
            expand_ranges, project_mmcif_object,
        )
        from openfold.data import mmcif_parsing  # noqa: E402

        print("\nconfirming flagged chains with the adapter's own check", flush=True)

        for n, chain in enumerate(flagged, 1):
            pdb, _, auth = chain.partition("_")
            entry = index[chain]
            try:
                parsed = mmcif_parsing.parse(
                    file_id=pdb,
                    mmcif_string=(args.mmcif_dir / f"{pdb}.cif").read_text(),
                )
                obj = parsed.mmcif_object
                if obj is None:
                    raise ValueError("mmcif_parsing returned no object")

                obj = project_mmcif_object(
                    obj, auth, expand_ranges(entry["r"]), entity_id=entry["e"])

                want = entry.get("s") or ""
                got = hashlib.sha256(
                    obj.chain_to_seqres[auth].encode()).hexdigest()[:len(want)]

                if got != want:
                    confirmed.append(
                        {"chain": chain, "expected": want, "observed": got})
            except Exception as exc:                      # noqa: BLE001
                confirmed.append({"chain": chain, "error": str(exc)[:160]})

            if n % 200 == 0:
                print(f"  confirmed {n:,} / {len(flagged):,}  "
                      f"bad so far: {len(confirmed)}", flush=True)

        print(f"\nCONFIRMED BAD: {len(confirmed):,} of {len(flagged):,} flagged")
        for c in confirmed[:20]:
            print(f"  {c}")

    args.output.write_text(json.dumps({
        "shard": args.shard,
        "num_shards": args.num_shards,
        "entries_scanned": len(mine),
        "flagged_chains": sorted(flagged),
        "flagged_count": len(flagged),
        "confirmed_bad": confirmed,
        "confirmed_bad_count": len(confirmed),
        "unreadable_entries": unreadable,
    }, indent=2, sort_keys=True))

    print(f"\nreport written  : {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
