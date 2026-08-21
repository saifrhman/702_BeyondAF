"""``pdbclean`` -- the single entry point for the COMP702 PDBClean pipeline.

The CLI resolves configuration, pins a snapshot, freezes a run identity, writes
provenance before expensive work begins, and then drives the *existing* stage
implementations.  It contains no scientific computation.

The UI (``pdbclean ui``) imports the same resolution and planning functions, so
a UI-configured run and a CLI-configured run produce the same
``resolved_run.yaml`` and execute the same backend.

Subcommands
-----------

``snapshots``       list the PDB snapshots available in the archive
``config``          show the fully resolved configuration and where each value came from
``plan``            show the stage plan for a configuration, including reuse decisions
``run``             create a run, write provenance, execute or submit the outstanding stages
``status``          list recorded runs, or show one run in detail
``stage-command``   print the argv a single stage would be executed with
``duplicates``      query the Duplicate Explorer from the terminal
``ui``              serve the web UI over the same backend
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
import sys
from pathlib import Path
from typing import Any, Sequence

from pdbclean import pipeline as pipeline_module
from pdbclean.duplicates import (
    DuplicateExplorer,
    DuplicateFilters,
    DuplicateQueryError,
    DuplicateSource,
)
from pdbclean.pipeline import (
    ACTION_BLOCKED,
    ACTION_REUSE,
    ACTION_RUN,
    PipelineError,
    PipelinePaths,
    build_executor,
    plan_pipeline,
    record_plan_in_provenance,
)
from pdbclean.run_provenance import RunProvenance, list_runs
from pdbclean.runconfig import (
    ResolvedRunConfig,
    RunConfigError,
    resolve_run_config,
)
from pdbclean.snapshot_selection import (
    SnapshotSelectionError,
    format_snapshot_id,
    interpret_menu_response,
    list_available_snapshots,
    render_snapshot_menu,
    resolve_snapshot_for_run,
    snapshot_provenance,
)
from pdbclean.stage_registry import STAGES_BY_ID, stage_catalogue


DEFAULT_PROFILE = "config/pdbclean/profiles/comp702_frozen_20260101.yaml"


def repository_root() -> Path:
    """Return the repository root.

    Derived from this module's location so a clone anywhere works, and
    overridable with ``PDBCLEAN_REPO_ROOT``.
    """

    override = os.environ.get("PDBCLEAN_REPO_ROOT")

    if override:
        return Path(override).resolve()

    return Path(__file__).resolve().parents[2]


# ----------------------------------------------------------------------
# Shared argument handling
# ----------------------------------------------------------------------


def add_config_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--config",
        default=None,
        help=(
            "Configuration file layered over the built-in validated defaults. "
            f"Use {DEFAULT_PROFILE} to reproduce the frozen COMP702 run."
        ),
    )
    parser.add_argument(
        "--set",
        dest="overrides",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help=(
            "Override one resolved value, e.g. "
            "--set duplicate_search.near_duplicate_threshold_angstrom=0.010. "
            "Repeatable. Overrides are recorded in provenance."
        ),
    )
    parser.add_argument(
        "--snapshot",
        default=None,
        help=(
            "Snapshot to use, as YYYY-MM-DD or YYYYMMDD. Omit to use the "
            "configured snapshot, or the latest complete snapshot."
        ),
    )
    parser.add_argument(
        "--repo-root",
        default=None,
        help="Repository root (defaults to the installed source tree).",
    )


def resolve_from_args(args: argparse.Namespace) -> ResolvedRunConfig:
    return resolve_run_config(
        config_path=args.config,
        overrides=args.overrides,
        override_origin="cli",
    )


def pin_snapshot(
    resolved: ResolvedRunConfig,
    args: argparse.Namespace,
    *,
    offline_ok: bool = True,
) -> tuple[ResolvedRunConfig, dict[str, Any]]:
    """Pin the snapshot, falling back to the configured identity when offline.

    Archive access is a network call.  When a concrete snapshot is already
    configured and the archive cannot be reached, the configured identity is
    used and the provenance records that it was not re-verified -- an
    inspection command must not fail merely because a login node is offline.
    A run, by contrast, always demands verification.
    """

    requested = args.snapshot

    try:
        pinned, snapshot = resolve_snapshot_for_run(
            resolved,
            requested=requested,
        )
        return pinned, snapshot_provenance(snapshot)
    except (SnapshotSelectionError, RunConfigError) as exc:
        configured = requested or resolved.get("snapshot.snapshot_id")

        if not offline_ok or not configured:
            raise

        from pdbclean.runconfig import with_resolved_snapshot
        from pdbclean.snapshot_selection import normalise_snapshot_id

        snapshot_id = normalise_snapshot_id(str(configured))

        pinned = with_resolved_snapshot(
            resolved,
            snapshot_id=snapshot_id,
            selection_mode="configured_unverified",
        )

        return pinned, {
            "snapshot_id": snapshot_id,
            "display": format_snapshot_id(snapshot_id),
            "selection_mode": "configured_unverified",
            "verified": False,
            "verification_error": str(exc),
        }


# ----------------------------------------------------------------------
# Rendering helpers
# ----------------------------------------------------------------------


def print_resolved_configuration(
    resolved: ResolvedRunConfig,
    *,
    stream=sys.stdout,
) -> None:
    print("Resolved scientific configuration", file=stream)
    print("=" * 78, file=stream)

    width = max(len(row["label"]) for row in resolved.scientific_summary())

    for row in resolved.scientific_summary():
        value = row["value"]
        rendered = "(unset)" if value is None else str(value)
        origin = row["source"]

        marker = "" if origin.startswith("builtin") else f"   [{origin}]"

        print(
            f"  {row['label']:<{width}}  {rendered}{marker}",
            file=stream,
        )

    print("", file=stream)
    print(f"  Resolved config SHA256   {resolved.sha256}", file=stream)
    print(f"  Scientific SHA256        {resolved.scientific_sha256}", file=stream)
    print(f"  Defaults version         {resolved.get('defaults_version')}", file=stream)

    if resolved.config_path:
        print(f"  Configuration file       {resolved.config_path}", file=stream)

    if resolved.override_items:
        print(
            f"  Explicit overrides       {', '.join(resolved.override_items)}",
            file=stream,
        )


STATUS_LABELS = {
    pipeline_module.PENDING: "pending",
    pipeline_module.BLOCKED: "blocked",
    pipeline_module.RUNNING: "running",
    pipeline_module.EXECUTION_COMPLETE: "exec complete",
    pipeline_module.VALIDATING: "validating",
    pipeline_module.PARTIAL: "partial",
    pipeline_module.VALIDATION_FAIL: "VALIDATION FAIL",
    pipeline_module.VALIDATION_PASS: "validation pass",
    pipeline_module.COMPLETE: "complete",
    pipeline_module.NOT_APPLICABLE: "n/a",
}


def print_plan(plan, *, stream=sys.stdout, verbose: bool = False) -> None:
    print(
        f"Pipeline plan   snapshot={plan.resolved.get('snapshot.snapshot_id')} "
        f"protocol={plan.resolved.get('release.protocol_version')}",
        file=stream,
    )
    print(f"Config SHA256   {plan.resolved.sha256}", file=stream)
    print(f"Science SHA256  {plan.resolved.scientific_sha256}", file=stream)
    print("=" * 100, file=stream)
    print(
        f"{'canonical':<15}{'stage':<26}{'layer':<9}{'status':<17}"
        f"{'validation':<17}{'action':<9}{'in':>11}{'out':>12}",
        file=stream,
    )

    for observation in plan.observations:
        stage = observation.stage

        print(
            f"{stage.canonical_stage:<15}{stage.stage_id:<26}{stage.layer:<9}"
            f"{STATUS_LABELS.get(observation.status, observation.status):<17}"
            f"{STATUS_LABELS.get(observation.validation, observation.validation):<17}"
            f"{observation.action:<9}"
            f"{_count(observation.input_count):>11}"
            f"{_count(observation.output_count):>12}",
            file=stream,
        )

        for issue in observation.incompatibilities:
            print(
                f"      ! {issue['summary_key']}: found {issue['observed']!r}, "
                f"configuration requires {issue['expected']!r}",
                file=stream,
            )

        if verbose:
            for message in observation.messages:
                print(f"      - {message}", file=stream)

    print("", file=stream)
    print(
        f"  reuse={len(plan.reusable)}  run={len(plan.to_run)}  "
        f"blocked={len(plan.blocked)}  complete={plan.complete}",
        file=stream,
    )


def _count(value: int | None) -> str:
    return "-" if value is None else f"{value:,}"


# ----------------------------------------------------------------------
# Stage command construction
# ----------------------------------------------------------------------


def stage_command(
    stage_id: str,
    resolved: ResolvedRunConfig,
    paths: PipelinePaths,
    *,
    protocol_config: str | None = None,
    python: str | None = None,
) -> list[str] | None:
    """Return the argv for one stage, or None when it has no command.

    Every scientific argument is derived from the resolved configuration, so
    the CLI, the UI and the Slurm wrappers all pass identical values.
    """

    interpreter = python or os.environ.get("PDBCLEAN_PYTHON") or sys.executable

    config_path = protocol_config or os.environ.get(
        "PDBCLEAN_PROTOCOL_CONFIG",
        str(paths.repo_root / "config/pdbclean/protocol_3_2_comp702_v1.yaml"),
    )

    policy_path = os.environ.get(
        "PDBCLEAN_REPRESENTATIVE_POLICY",
        str(
            paths.repo_root
            / "config/pdbclean/stage14_representative_policy_v1.yaml"
        ),
    )

    protocol_root = paths.output_root / paths.snapshot / paths.protocol
    threshold_mA = resolved.near_duplicate_threshold_mA

    expectations = resolved.get("expectations") or {}

    if stage_id == "redundancy_graph":
        full_bri = protocol_root / "full_bri_nn" / "finalized"

        command = [
            interpreter,
            str(paths.repo_root / "scripts/build_stage14_geometric_graph.py"),
            "--edges",
            str(full_bri / "candidate_near_duplicates.parquet"),
            "--m1-edges",
            str(full_bri / "m1_near_duplicates.parquet"),
            "--output-dir",
            str(protocol_root / "stage14_geometric_graph"),
            "--threshold-mA",
            str(threshold_mA),
            "--minimum-chain-length",
            str(
                resolved.get(
                    "representative_selection."
                    "minimum_deduplicated_chain_length"
                )
            ),
        ]

        # The Stage-6 summary lets a snapshot whose counts are not known in
        # advance derive its own m >= 2 population.
        command.extend(
            [
                "--length-buckets-summary",
                str(protocol_root / "length_buckets" / "global_summary.json"),
            ]
        )

        # Dataset-version gates are asserted when this snapshot's counts are
        # known.  A new snapshot has none, and must not be required to match
        # the 2026-01-01 figures.
        supplied = False

        for flag, key in (
            ("--expected-edges", "edge_count"),
            ("--expected-m1-edges", "m1_edge_count"),
            ("--expected-mge2-nodes", "mge2_node_count"),
        ):
            value = expectations.get(key)

            if value is not None:
                command.extend([flag, str(value)])
                supplied = True

        if not supplied:
            command.append("--no-expectation-gate")

        return command

    if stage_id == "representative_selection":
        command = [
            interpreter,
            str(paths.repo_root / "scripts/select_stage14_representatives.py"),
            "--graph-dir",
            str(protocol_root / "stage14_geometric_graph"),
            "--edges",
            str(
                protocol_root
                / "full_bri_nn/finalized/candidate_near_duplicates.parquet"
            ),
            "--accepted",
            str(protocol_root / "quality/merged/accepted.parquet"),
            "--metadata",
            str(
                protocol_root
                / "downstream_metadata/finalized/entry_metadata.parquet"
            ),
            "--config",
            policy_path,
            "--output-dir",
            str(protocol_root / "stage14_representative_selection_v1"),
            "--threshold-mA",
            str(threshold_mA),
        ]

        canonical = expectations.get("canonical_input_chain_count")

        if canonical is not None:
            command.extend(
                ["--expected-canonical-input-chains", str(canonical)]
            )

        return command

    if stage_id == "gold_release":
        command = [
            interpreter,
            str(paths.repo_root / "scripts/build_stage14_final_release.py"),
            "--protocol-root",
            str(protocol_root),
            "--policy-config",
            policy_path,
            "--output-dir",
            str(paths.release_root / paths.release),
            "--threshold-mA",
            str(threshold_mA),
        ]

        retained = expectations.get("retained_chain_count")
        removed = expectations.get("removed_chain_count")

        if retained is not None and removed is not None:
            command.extend(
                [
                    "--expected-retained-chains",
                    str(retained),
                    "--expected-removed-chains",
                    str(removed),
                ]
            )
        else:
            command.append("--no-expectation-gate")

        return command

    if stage_id == "candidate_filtering":
        return [
            interpreter,
            "-m",
            "pdbclean.brain_prefilter_production",
            "--config",
            config_path,
        ]

    if stage_id == "complete_bri_nn":
        return [
            interpreter,
            "-m",
            "pdbclean.full_bri_nn_production",
            "--config",
            config_path,
        ]

    if stage_id == "duplicate_classification":
        return [
            interpreter,
            "-m",
            "pdbclean.duplicate_classification_production",
            "--config",
            config_path,
        ]

    if stage_id == "length_buckets":
        return [
            interpreter,
            "-m",
            "pdbclean.length_buckets_cli",
            "--config",
            config_path,
        ]

    if stage_id == "brain":
        return [
            interpreter,
            "-m",
            "pdbclean.brain_finalize_cli",
            "--config",
            config_path,
        ]

    if stage_id == "structural_cleaning":
        return [
            "bash",
            str(paths.repo_root / "scripts/pdbclean/submit_quality_pipeline.sh"),
            config_path,
            str(
                paths.output_root
                / paths.snapshot
                / "bronze/source_manifest.parquet"
            ),
        ]

    if stage_id == "bronze_source_manifest":
        return [
            interpreter,
            str(paths.repo_root / "scripts/pdbclean/create_manifest.py"),
            "--config",
            config_path,
        ]

    return None


# ----------------------------------------------------------------------
# Subcommands
# ----------------------------------------------------------------------


def cmd_snapshots(args: argparse.Namespace) -> int:
    resolved = resolve_from_args(args)

    bucket_url = resolved.get("snapshot.bucket_url")

    try:
        choices = list_available_snapshots(
            bucket_url=bucket_url,
            limit=args.limit,
        )
    except SnapshotSelectionError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if args.json:
        print(
            json.dumps(
                [
                    {
                        "index": choice.index,
                        "snapshot_id": choice.snapshot_id,
                        "display": choice.display,
                        "is_latest": choice.is_latest,
                    }
                    for choice in choices
                ],
                indent=2,
            )
        )
        return 0

    print(render_snapshot_menu(choices))
    print("")
    print(f"Archive: {bucket_url}")
    print(f"Discovered: {len(choices)} dated snapshots")

    return 0


def cmd_config(args: argparse.Namespace) -> int:
    resolved = resolve_from_args(args)

    if args.snapshot or resolved.get("snapshot.snapshot_id"):
        resolved, _ = pin_snapshot(resolved, args)

    if args.json:
        print(json.dumps(resolved.to_dict(), indent=2, sort_keys=True))
        return 0

    if args.yaml:
        print(resolved.to_yaml(), end="")
        return 0

    print_resolved_configuration(resolved)

    if args.sources:
        print("")
        print("Value provenance")
        print("=" * 78)

        for key in sorted(resolved.sources):
            print(f"  {key:<64} {resolved.sources[key]}")

    return 0


def cmd_plan(args: argparse.Namespace) -> int:
    resolved = resolve_from_args(args)
    resolved, _ = pin_snapshot(resolved, args)

    repo_root = Path(args.repo_root) if args.repo_root else repository_root()

    plan = plan_pipeline(
        resolved,
        repo_root=repo_root,
        compute_checksums=args.checksums,
    )

    if args.json:
        print(json.dumps(plan.to_dict(), indent=2, sort_keys=True))
        return 0

    print_plan(plan, verbose=args.verbose)

    return 0


def cmd_stage_command(args: argparse.Namespace) -> int:
    resolved = resolve_from_args(args)
    resolved, _ = pin_snapshot(resolved, args)

    repo_root = Path(args.repo_root) if args.repo_root else repository_root()
    paths = PipelinePaths.from_config(resolved, repo_root=repo_root)

    if args.stage not in STAGES_BY_ID:
        print(f"error: unknown stage {args.stage!r}", file=sys.stderr)
        return 2

    try:
        command = stage_command(
            args.stage,
            resolved,
            paths,
            protocol_config=args.protocol_config,
        )
    except PipelineError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if command is None:
        print(
            f"error: stage {args.stage!r} has no executable command",
            file=sys.stderr,
        )
        return 2

    if args.shell:
        print(" ".join(shlex.quote(part) for part in command))
    else:
        for part in command:
            print(part)

    return 0


def cmd_run(args: argparse.Namespace) -> int:
    resolved = resolve_from_args(args)

    repo_root = Path(args.repo_root) if args.repo_root else repository_root()

    # ---- snapshot selection -------------------------------------------
    if args.interactive:
        bucket_url = resolved.get("snapshot.bucket_url")

        try:
            choices = list_available_snapshots(
                bucket_url=bucket_url,
                limit=args.limit,
            )
        except SnapshotSelectionError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2

        print(render_snapshot_menu(choices))
        print("")

        response = input("Select snapshot [1]: ")

        try:
            selected = interpret_menu_response(response, choices)
        except SnapshotSelectionError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2

        args.snapshot = selected

    try:
        resolved, snapshot = pin_snapshot(
            resolved,
            args,
            offline_ok=False,
        )
    except (SnapshotSelectionError, RunConfigError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    paths = PipelinePaths.from_config(resolved, repo_root=repo_root)

    # ---- show the resolved configuration before anything expensive ----
    print_resolved_configuration(resolved)
    print("")
    print(f"  Snapshot                 {snapshot.get('display')}")
    print(f"  Release name             {paths.release}")
    print(f"  Output root              {paths.output_root}")
    print(f"  Release root             {paths.release_root}")
    print(f"  Executor                 {args.executor}")
    print("")

    plan = plan_pipeline(resolved, repo_root=repo_root)

    print_plan(plan)
    print("")

    if args.plan_only:
        return 0

    if not args.yes:
        response = input("Start this run? [y/N]: ").strip().lower()

        if response not in {"y", "yes"}:
            print("Aborted; nothing was executed and no run was created.")
            return 1

    # ---- freeze run identity and write provenance BEFORE work ---------
    provenance = RunProvenance.create(
        resolved=resolved,
        run_root=paths.run_root,
        repo_root=repo_root,
        snapshot=snapshot,
        invocation={
            "argv": sys.argv[1:],
            "executor": args.executor,
            "origin": "cli",
        },
    )

    print(f"Run ID          {provenance.run_id}")
    print(f"Run directory   {provenance.run_dir}")
    print("")

    record_plan_in_provenance(plan, provenance)

    executor = build_executor(args.executor)

    provenance.set_status("running", executor=executor.name)

    failed = False

    for observation in plan.observations:
        stage = observation.stage

        if observation.action == ACTION_REUSE:
            print(
                f"  reuse    {stage.stage_id:<26} "
                f"(validated output already present)"
            )
            continue

        if observation.action == pipeline_module.ACTION_NOT_APPLICABLE:
            continue

        if observation.action == ACTION_BLOCKED:
            print(
                f"  BLOCKED  {stage.stage_id:<26} "
                f"upstream validation has not passed"
            )
            failed = True
            break

        try:
            command = stage_command(stage.stage_id, resolved, paths)
        except PipelineError as exc:
            print(f"  ERROR    {stage.stage_id:<26} {exc}")
            provenance.update_stage(
                stage.stage_id,
                status=pipeline_module.VALIDATION_FAIL,
                messages=[str(exc)],
            )
            failed = True
            break

        if command is None:
            print(
                f"  skip     {stage.stage_id:<26} "
                f"no executable command registered"
            )
            continue

        rendered = " ".join(shlex.quote(part) for part in command)

        if executor.name == "dry-run":
            print(f"  would run {stage.stage_id:<25} {rendered}")
            continue

        print(f"  run      {stage.stage_id:<26} {rendered}")

        # Host-specific facts belong to the execution, not to the run's
        # canonical configuration.
        provenance.record_runtime(resolved, stage_id=stage.stage_id)

        provenance.update_stage(
            stage.stage_id,
            status=pipeline_module.RUNNING,
            attempts=(observation_attempts(provenance, stage.stage_id) + 1),
        )

        result = executor.run(
            stage.stage_id,
            command,
            cwd=paths.repo_root,
        )

        provenance.append_event("stage_execution", **result.to_dict())

        if result.slurm_job_id:
            stage_record = provenance.stage(stage.stage_id)

            if stage_record is not None:
                stage_record.slurm_job_ids.append(result.slurm_job_id)

            print(f"           submitted as Slurm job {result.slurm_job_id}")

        if result.returncode not in (None, 0):
            provenance.update_stage(
                stage.stage_id,
                status=pipeline_module.VALIDATION_FAIL,
                validation=pipeline_module.VALIDATION_FAIL,
            )
            provenance.record_validation(
                f"{stage.stage_id}:execution",
                "FAIL",
                returncode=result.returncode,
            )
            print(
                f"  FAILED   {stage.stage_id:<26} "
                f"exit {result.returncode}; downstream stages will not start"
            )
            failed = True
            break

        provenance.update_stage(
            stage.stage_id,
            status=pipeline_module.EXECUTION_COMPLETE,
        )

    provenance.set_status("failed" if failed else "completed")
    provenance.flush()

    print("")
    print(f"Provenance      {provenance.record_path}")

    return 1 if failed else 0


def observation_attempts(provenance: RunProvenance, stage_id: str) -> int:
    stage = provenance.stage(stage_id)

    return 0 if stage is None else stage.attempts


def cmd_status(args: argparse.Namespace) -> int:
    resolved = resolve_from_args(args)
    repo_root = Path(args.repo_root) if args.repo_root else repository_root()

    run_root = resolved.get("storage.run_root") or "outputs/runs"
    run_path = Path(run_root)

    if not run_path.is_absolute():
        run_path = repo_root / run_path

    if args.run_id:
        directory = run_path / args.run_id

        if not directory.is_dir():
            print(f"error: no such run: {directory}", file=sys.stderr)
            return 2

        record = RunProvenance.load(directory).record

        if args.json:
            print(json.dumps(record, indent=2, sort_keys=True))
            return 0

        print(f"Run          {record['run_id']}")
        print(f"Status       {record.get('status')}")
        print(f"Created      {record.get('created_at')}")
        print(f"Snapshot     {(record.get('snapshot') or {}).get('display')}")
        print(f"Config sha   {record.get('resolved_config_sha256')}")
        print(f"Science sha  {record.get('scientific_config_sha256')}")
        print(f"Git branch   {(record.get('git') or {}).get('branch')}")
        print(f"Git commit   {(record.get('git') or {}).get('commit')}")
        print(f"Tree dirty   {(record.get('git') or {}).get('working_tree_dirty')}")
        print("")
        print(f"{'stage':<28}{'status':<18}{'validation':<18}{'out':>12}")

        for stage in record.get("stages", []):
            print(
                f"{stage['stage_id']:<28}{stage.get('status', ''):<18}"
                f"{stage.get('validation', ''):<18}"
                f"{_count(stage.get('output_count')):>12}"
            )

        return 0

    runs = list_runs(run_path)

    if args.json:
        print(json.dumps(runs, indent=2, sort_keys=True))
        return 0

    if not runs:
        print(f"No runs recorded under {run_path}")
        return 0

    print(f"{'run id':<34}{'created':<22}{'status':<12}{'snapshot':<12}")

    for run in runs:
        print(
            f"{run['run_id']:<34}{str(run.get('created_at')):<22}"
            f"{str(run.get('status')):<12}{str(run.get('snapshot_id')):<12}"
        )

    return 0


def cmd_duplicates(args: argparse.Namespace) -> int:
    resolved = resolve_from_args(args)
    resolved, _ = pin_snapshot(resolved, args)

    repo_root = Path(args.repo_root) if args.repo_root else repository_root()
    paths = PipelinePaths.from_config(resolved, repo_root=repo_root)

    source = DuplicateSource(
        protocol_root=paths.output_root / paths.snapshot / paths.protocol,
        release_root=paths.release_root / paths.release,
    )

    try:
        explorer = DuplicateExplorer(source)
    except DuplicateQueryError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if args.summary:
        print(json.dumps(explorer.summary(), indent=2, sort_keys=True))
        return 0

    filters = DuplicateFilters(
        pdb_id=args.pdb_id,
        chain=args.chain,
        exact_only=args.exact_only,
        nonzero_near_only=args.nonzero_near_only,
        min_length=args.min_length,
        max_length=args.max_length,
        min_distance_mA=args.min_distance,
        max_distance_mA=args.max_distance,
        relationship=args.relationship,
        offset=args.offset,
        limit=args.limit,
    )

    try:
        result = explorer.query(filters)
    except DuplicateQueryError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
        return 0

    print(
        f"{result['matched']:,} matching pairs; showing "
        f"{len(result['rows'])} from offset {result['offset']}"
    )
    print("")
    print(
        f"{'A':<12}{'B':<12}{'len':>7}{'d(mA)':>8}{'d(A)':>10}  "
        f"{'class':<22}{'relationship':<14}{'representative':<14}"
    )

    for row in result["rows"]:
        print(
            f"{row['pdb_id_a'] + ':' + row['chain_a']:<12}"
            f"{row['pdb_id_b'] + ':' + row['chain_b']:<12}"
            f"{row['chain_length']:>7}{row['d_bri_mA']:>8}"
            f"{row['d_bri_angstrom']:>10.4f}  "
            f"{row['classification']:<22}{row['relationship']:<14}"
            f"{str(row['representative'] or '-'):<14}"
        )

    return 0


def cmd_stages(args: argparse.Namespace) -> int:
    catalogue = stage_catalogue()

    if args.json:
        print(json.dumps(catalogue, indent=2))
        return 0

    for stage in catalogue:
        print(
            f"{stage['canonical_stage']:<15}{stage['title']}  "
            f"[{stage['layer']}]  ({stage['stage_id']})"
        )
        print(f"    {stage['purpose']}")

        if stage["validation"]:
            print(f"    gate: {stage['validation']}")

        print("")

    return 0


def cmd_ui(args: argparse.Namespace) -> int:
    from pdbclean.ui.server import serve

    repo_root = Path(args.repo_root) if args.repo_root else repository_root()

    return serve(
        host=args.host,
        port=args.port,
        repo_root=repo_root,
        config_path=args.config,
        overrides=args.overrides,
        open_browser=not args.no_browser,
    )


# ----------------------------------------------------------------------
# Parser
# ----------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pdbclean",
        description=(
            "COMP702 PDBClean: geometric redundancy detection and removal "
            "over a PDB snapshot."
        ),
    )

    sub = parser.add_subparsers(dest="command", required=True)

    # snapshots ---------------------------------------------------------
    p = sub.add_parser("snapshots", help="List available PDB snapshots")
    add_config_arguments(p)
    p.add_argument("--limit", type=int, default=20)
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_snapshots)

    # config ------------------------------------------------------------
    p = sub.add_parser("config", help="Show the resolved configuration")
    add_config_arguments(p)
    p.add_argument("--json", action="store_true")
    p.add_argument("--yaml", action="store_true")
    p.add_argument(
        "--sources",
        action="store_true",
        help="Show which layer supplied every resolved value",
    )
    p.set_defaults(func=cmd_config)

    # plan --------------------------------------------------------------
    p = sub.add_parser("plan", help="Show the stage plan")
    add_config_arguments(p)
    p.add_argument("--json", action="store_true")
    p.add_argument("--verbose", "-v", action="store_true")
    p.add_argument(
        "--checksums",
        action="store_true",
        help="Compute output checksums (slower)",
    )
    p.set_defaults(func=cmd_plan)

    # run ---------------------------------------------------------------
    p = sub.add_parser("run", help="Create a run and execute outstanding stages")
    add_config_arguments(p)
    p.add_argument(
        "--executor",
        default="dry-run",
        choices=["dry-run", "local", "slurm"],
        help=(
            "How stages are realised. dry-run prints commands and is the "
            "default; slurm submits batch scripts with sbatch."
        ),
    )
    p.add_argument(
        "--interactive",
        action="store_true",
        help="Choose the snapshot from a menu before starting",
    )
    p.add_argument("--limit", type=int, default=20)
    p.add_argument("--plan-only", action="store_true")
    p.add_argument(
        "--yes",
        "-y",
        action="store_true",
        help="Do not prompt for confirmation",
    )
    p.set_defaults(func=cmd_run)

    # status ------------------------------------------------------------
    p = sub.add_parser("status", help="Show recorded runs")
    add_config_arguments(p)
    p.add_argument("run_id", nargs="?", default=None)
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_status)

    # stage-command -----------------------------------------------------
    p = sub.add_parser(
        "stage-command",
        help="Print the argv for one stage",
    )
    add_config_arguments(p)
    p.add_argument("--stage", required=True)
    p.add_argument("--protocol-config", default=None)
    p.add_argument("--shell", action="store_true")
    p.set_defaults(func=cmd_stage_command)

    # stages ------------------------------------------------------------
    p = sub.add_parser("stages", help="Describe the pipeline stages")
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_stages)

    # duplicates --------------------------------------------------------
    p = sub.add_parser("duplicates", help="Query detected duplicate pairs")
    add_config_arguments(p)
    p.add_argument("--pdb-id", default=None)
    p.add_argument("--chain", default=None)
    p.add_argument("--exact-only", action="store_true")
    p.add_argument("--nonzero-near-only", action="store_true")
    p.add_argument("--min-length", type=int, default=None)
    p.add_argument("--max-length", type=int, default=None)
    p.add_argument("--min-distance", type=int, default=None, metavar="MA")
    p.add_argument("--max-distance", type=int, default=None, metavar="MA")
    p.add_argument(
        "--relationship",
        default=None,
        choices=["removed", "retained", "unaffected"],
    )
    p.add_argument("--offset", type=int, default=0)
    p.add_argument("--limit", type=int, default=25)
    p.add_argument("--summary", action="store_true")
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_duplicates)

    # ui ----------------------------------------------------------------
    p = sub.add_parser("ui", help="Serve the web UI")
    add_config_arguments(p)
    p.add_argument("--host", default="127.0.0.1")
    p.add_argument("--port", type=int, default=8765)
    p.add_argument("--no-browser", action="store_true")
    p.set_defaults(func=cmd_ui)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        return args.func(args)
    except (RunConfigError, PipelineError, SnapshotSelectionError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print("\ninterrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
