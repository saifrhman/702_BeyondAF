"""HTTP backend for the PDBClean UI.

Implemented on the Python standard library so the UI adds no dependency to the
pipeline environment.  It is a local research tool: it binds to loopback by
default and serves only whitelisted directories.

Every endpoint delegates to the same modules the CLI uses.  A UI-configured run
and a CLI-configured run resolve to the same ``resolved_run.yaml`` and execute
the same stage commands.
"""

from __future__ import annotations

import json
import mimetypes
import posixpath
import threading
import traceback
import webbrowser
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, quote, unquote, urlparse

from pdbclean import pipeline as pipeline_module
from pdbclean.defaults import DEFAULTS_VERSION, validated_defaults
from pdbclean.duplicates import (
    DuplicateExplorer,
    DuplicateFilters,
    DuplicateQueryError,
    DuplicateSource,
)
from pdbclean.pipeline import (
    PipelineError,
    PipelinePaths,
    gold_release_summary,
    plan_pipeline,
    record_plan_in_provenance,
)
from pdbclean import artefacts as artefact_module
from pdbclean.run_inspection import (
    duplicate_navigation,
    run_timeline,
    stage_detail,
)
from pdbclean.molstar_service import (
    MolstarSceneService,
    MolstarServiceError,
    PairRequest,
    locator_for_run,
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
    list_available_snapshots,
)
from pdbclean.defaults import (
    IMPLEMENTED_REPRESENTATION_PRECISION_ANGSTROM,
    precision_is_implemented,
    representation_unit_label,
)
from pdbclean.snapshot_store import SnapshotStoreLayout, snapshot_status
from pdbclean.stage_registry import (
    LAYERS,
    canonical_catalogue,
    canonical_for_producer,
    canonical_timeline,
    stage_catalogue,
)


STATIC_DIR = Path(__file__).resolve().parent / "static"

#: Directories the UI may serve structure files from, relative to the repo.
STRUCTURE_ROOTS = (
    "reports/molstar_exact_duplicate_examples",
)

#: File suffixes the structure route will serve.
STRUCTURE_SUFFIXES = frozenset(
    {".cif", ".mvsj", ".json", ".bcif", ".pdb", ".mmcif"}
)


#: Largest artefact the server will stream to a browser. Above this the UI
#: reports the path so the file can be copied directly instead.
MAX_DOWNLOAD_BYTES = 512 * 1024 * 1024


class UIError(RuntimeError):
    """Raised for a request the UI cannot satisfy."""


class ArtefactForbidden(UIError):
    """Raised when a requested path lies outside the allowed roots."""


class UIState:
    """Process-wide state shared by request handlers."""

    def __init__(
        self,
        *,
        repo_root: Path,
        config_path: str | None,
        overrides: list[str],
    ) -> None:
        self.repo_root = repo_root
        self.default_config_path = config_path
        self.default_overrides = list(overrides or [])
        self._lock = threading.Lock()
        self._explorers: dict[str, DuplicateExplorer] = {}

    # -- configuration --------------------------------------------------

    def resolve(
        self,
        *,
        config_path: str | None = None,
        overrides: Any = None,
        snapshot: str | None = None,
    ) -> ResolvedRunConfig:
        resolved = resolve_run_config(
            config_path=(
                config_path
                if config_path is not None
                else self.default_config_path
            ),
            overrides=(
                overrides
                if overrides is not None
                else self.default_overrides
            ),
            override_origin="ui",
        )

        if snapshot:
            from pdbclean.runconfig import with_resolved_snapshot
            from pdbclean.snapshot_selection import normalise_snapshot_id

            resolved = with_resolved_snapshot(
                resolved,
                snapshot_id=normalise_snapshot_id(snapshot),
                selection_mode="ui_selection",
            )

        return resolved

    def scene_cache_root(self, resolved: ResolvedRunConfig) -> Path:
        """Disposable Mol* scene cache. Never a scientific artefact."""

        configured = resolved.get("storage.hot_cache_root") or (
            "outputs/snapshot_cache"
        )

        root = Path(configured)

        if not root.is_absolute():
            root = self.repo_root / root

        return root / "molstar_scenes"

    def molstar(
        self,
        resolved: ResolvedRunConfig,
        *,
        snapshot_id: str | None = None,
    ) -> MolstarSceneService:
        """Scene service scoped to one run's snapshot."""

        snapshot = snapshot_id or resolved.get("snapshot.snapshot_id")

        def _root(dotted: str, fallback: str) -> Path:
            value = resolved.get(dotted) or fallback
            candidate = Path(value)

            return (
                candidate
                if candidate.is_absolute()
                else self.repo_root / candidate
            )

        return MolstarSceneService(
            locator=locator_for_run(
                repo_root=self.repo_root,
                snapshot_id=str(snapshot) if snapshot else None,
                hot_cache_root=_root(
                    "storage.hot_cache_root", "outputs/snapshot_cache"
                ),
                durable_root=_root(
                    "storage.durable_snapshot_root", "outputs/snapshot_store"
                ),
                # The run's own Bronze manifest supplies the immutable
                # source-object identity for this snapshot.
                output_root=_root("storage.output_root", "outputs/pdbclean"),
                bucket_url=resolved.get("snapshot.bucket_url"),
                protocol=resolved.get("release.protocol_version"),
            ),
            cache_root=self.scene_cache_root(resolved),
        )

    def paths(self, resolved: ResolvedRunConfig) -> PipelinePaths:
        return PipelinePaths.from_config(resolved, repo_root=self.repo_root)

    # -- duplicate explorer ---------------------------------------------

    def explorer(self, resolved: ResolvedRunConfig) -> DuplicateExplorer:
        paths = self.paths(resolved)

        key = f"{paths.snapshot}:{paths.protocol}"

        with self._lock:
            existing = self._explorers.get(key)

            if existing is not None:
                return existing

            explorer = DuplicateExplorer(
                DuplicateSource(
                    protocol_root=(
                        paths.output_root / paths.snapshot / paths.protocol
                    ),
                    release_root=paths.release_root / paths.release,
                )
            )

            self._explorers[key] = explorer

            return explorer


def _json_bytes(payload: Any) -> bytes:
    return json.dumps(payload, indent=2, sort_keys=True).encode("utf-8")


class Handler(BaseHTTPRequestHandler):
    server_version = "PDBCleanUI/1.0"
    state: UIState

    # Quieter, structured request logging.
    def log_message(self, fmt: str, *args: Any) -> None:  # noqa: A003
        print(f"[ui] {self.address_string()} {fmt % args}")

    # -- plumbing -------------------------------------------------------

    def _send(
        self,
        status: HTTPStatus,
        body: bytes,
        content_type: str = "application/json; charset=utf-8",
        *,
        extra_headers: dict[str, str] | None = None,
    ) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Cache-Control", "no-store")

        for key, value in (extra_headers or {}).items():
            self.send_header(key, value)

        self.end_headers()

        if self.command != "HEAD":
            self.wfile.write(body)

    def _json(self, payload: Any, status: HTTPStatus = HTTPStatus.OK) -> None:
        self._send(status, _json_bytes(payload))

    def _error(
        self,
        message: str,
        status: HTTPStatus = HTTPStatus.BAD_REQUEST,
        **extra: Any,
    ) -> None:
        self._json({"error": message, **extra}, status)

    def _read_body(self) -> dict[str, Any]:
        length = int(self.headers.get("Content-Length") or 0)

        if length <= 0:
            return {}

        raw = self.rfile.read(length)

        try:
            payload = json.loads(raw.decode("utf-8"))
        except (UnicodeDecodeError, ValueError) as exc:
            raise UIError(f"Request body is not valid JSON: {exc}") from exc

        if not isinstance(payload, dict):
            raise UIError("Request body must be a JSON object")

        return payload

    # -- routing --------------------------------------------------------

    def do_GET(self) -> None:  # noqa: N802
        parsed = urlparse(self.path)
        route = parsed.path
        query = parse_qs(parsed.query)

        try:
            if route.startswith("/api/"):
                self._route_api_get(route, query)
                return

            if route.startswith("/structures/"):
                self._serve_structure(route)
                return

            self._serve_static(route)
        except UIError as exc:
            self._error(str(exc))
        except (RunConfigError, PipelineError, SnapshotSelectionError) as exc:
            self._error(str(exc))
        except DuplicateQueryError as exc:
            self._error(str(exc), HTTPStatus.NOT_FOUND)
        except Exception as exc:  # pragma: no cover - defensive
            traceback.print_exc()
            self._error(
                f"{type(exc).__name__}: {exc}",
                HTTPStatus.INTERNAL_SERVER_ERROR,
            )

    def do_HEAD(self) -> None:  # noqa: N802
        self.do_GET()

    def do_POST(self) -> None:  # noqa: N802
        parsed = urlparse(self.path)

        try:
            body = self._read_body()

            if parsed.path == "/api/config/resolve":
                self._api_resolve(body)
                return

            if parsed.path == "/api/plan":
                self._api_plan(body)
                return

            if parsed.path == "/api/run":
                self._api_run(body)
                return

            self._error("Unknown endpoint", HTTPStatus.NOT_FOUND)
        except UIError as exc:
            self._error(str(exc))
        except (RunConfigError, PipelineError, SnapshotSelectionError) as exc:
            self._error(str(exc))
        except Exception as exc:  # pragma: no cover - defensive
            traceback.print_exc()
            self._error(
                f"{type(exc).__name__}: {exc}",
                HTTPStatus.INTERNAL_SERVER_ERROR,
            )

    # -- GET endpoints --------------------------------------------------

    def _route_api_get(self, route: str, query: dict[str, list[str]]) -> None:
        if route == "/api/bootstrap":
            self._api_bootstrap()
            return

        if route == "/api/snapshots":
            self._api_snapshots(query)
            return

        if route == "/api/duplicates":
            self._api_duplicates(query)
            return

        if route == "/api/release":
            self._api_release(query)
            return

        if route == "/api/runs":
            self._api_runs(query)
            return

        if route == "/api/artefact":
            self._api_artefact(query)
            return

        if route == "/api/artefact/table":
            self._api_artefact_table(query)
            return

        if route == "/api/artefact/download":
            self._api_artefact_download(query)
            return

        if route == "/api/structure":
            self._api_structure(query)
            return

        if route == "/api/pair/availability":
            self._api_pair_availability(query)
            return

        if route == "/api/pair/scene":
            self._api_pair_scene(query)
            return

        if route.startswith("/api/runs/"):
            remainder = route[len("/api/runs/"):]
            parts = [part for part in remainder.split("/") if part]

            if len(parts) == 1:
                self._api_run_detail(parts[0])
                return

            if len(parts) == 2 and parts[1] == "timeline":
                self._api_run_timeline(parts[0])
                return

            if len(parts) == 3 and parts[1] == "stages":
                self._api_run_stage(parts[0], parts[2], query)
                return

            self._error("Unknown run endpoint", HTTPStatus.NOT_FOUND)
            return

        if route == "/api/scenes":
            self._api_scenes()
            return

        self._error("Unknown endpoint", HTTPStatus.NOT_FOUND)

    def _api_bootstrap(self) -> None:
        state = self.state

        profiles_dir = state.repo_root / "config" / "pdbclean" / "profiles"

        profiles = []

        if profiles_dir.is_dir():
            for path in sorted(profiles_dir.glob("*.yaml")):
                profiles.append(
                    {
                        "path": str(path.relative_to(state.repo_root)),
                        "name": path.stem,
                    }
                )

        self._json(
            {
                "repo_root": str(state.repo_root),
                "defaults_version": DEFAULTS_VERSION,
                "defaults": validated_defaults(),
                "default_config_path": state.default_config_path,
                "profiles": profiles,
                "stages": stage_catalogue(),
                "canonical_timeline": canonical_catalogue(),
                "layers": [
                    {"id": layer, "label": label} for layer, label in LAYERS
                ],
                "stage_states": list(pipeline_module.STATE_ORDER),
            }
        )

    def _api_snapshots(self, query: dict[str, list[str]]) -> None:
        resolved = self.state.resolve(
            config_path=_single(query, "config"),
        )

        limit = int(_single(query, "limit") or 20)

        try:
            choices = list_available_snapshots(
                bucket_url=resolved.get("snapshot.bucket_url"),
                limit=limit,
            )
        except SnapshotSelectionError as exc:
            self._json(
                {
                    "bucket_url": resolved.get("snapshot.bucket_url"),
                    "snapshots": [],
                    "error": str(exc),
                }
            )
            return

        self._json(
            {
                "bucket_url": resolved.get("snapshot.bucket_url"),
                "snapshots": [
                    {
                        "index": choice.index,
                        "snapshot_id": choice.snapshot_id,
                        "display": choice.display,
                        "is_latest": choice.is_latest,
                    }
                    for choice in choices
                ],
            }
        )

    def _api_duplicates(self, query: dict[str, list[str]]) -> None:
        resolved = self.state.resolve(
            config_path=_single(query, "config"),
            snapshot=_single(query, "snapshot"),
        )

        explorer = self.state.explorer(resolved)

        filters = DuplicateFilters(
            pdb_id=_single(query, "pdb_id"),
            chain=_single(query, "chain"),
            exact_only=_flag(query, "exact_only"),
            nonzero_near_only=_flag(query, "nonzero_near_only"),
            min_length=_int(query, "min_length"),
            max_length=_int(query, "max_length"),
            min_distance_mA=_int(query, "min_distance"),
            max_distance_mA=_int(query, "max_distance"),
            relationship=_single(query, "relationship"),
            offset=_int(query, "offset") or 0,
            limit=_int(query, "limit") or 50,
        )

        result = explorer.query(filters)
        result["summary"] = explorer.summary()
        result["scenes"] = _scene_index(self.state.repo_root)

        # Per-row Mol* availability, resolved against THIS run's snapshot.
        # Every row gets either a usable action or a specific reason -- never
        # a bare "no prepared scene".
        service = self.state.molstar(resolved)
        cache: dict[tuple[str, str], dict[str, Any]] = {}

        for row in result["rows"]:
            request = PairRequest(
                pdb_id_a=row["pdb_id_a"],
                chain_a=row["chain_a"],
                pdb_id_b=row["pdb_id_b"],
                chain_b=row["chain_b"],
                snapshot_id=str(resolved.get("snapshot.snapshot_id") or ""),
                model_id=row.get("model_a", 1),
                chain_length=row.get("chain_length"),
                d_bri_units=row.get("d_bri_mA"),
                classification=row.get("classification"),
                relationship=row.get("relationship"),
                representative=row.get("representative"),
            )

            key = (row["pdb_id_a"], row["pdb_id_b"])

            if key not in cache:
                availability = service.availability(request)
                cache[key] = {
                    "available": availability["available"],
                    "reason": availability["reason"],
                    "reason_text": availability["reason_text"],
                    "views": [
                        view["key"] for view in availability.get("views", [])
                    ],
                }

            row["molstar"] = dict(cache[key])
            row["molstar"]["query"] = {
                "pdb_id_a": request.pdb_id_a,
                "chain_a": request.chain_a,
                "pdb_id_b": request.pdb_id_b,
                "chain_b": request.chain_b,
                "snapshot": request.snapshot_id,
                "model_id": request.model_id,
                "d_bri": request.d_bri_units,
                "classification": request.classification,
                "relationship": request.relationship,
                "representative": request.representative or "",
            }

        # The authoritative tables behind these rows, openable in the viewer.
        result["source_tables"] = [
            {
                "path": path,
                "name": Path(path).name,
                "role": role,
            }
            for role, paths in (
                ("classification", explorer.source.classification_paths),
                ("near_duplicates", explorer.source.near_duplicate_paths),
            )
            for path in [str(candidate) for candidate in paths]
            if Path(path).is_file()
        ]

        mapping = explorer.source.representative_mapping_path

        if mapping is not None:
            result["source_tables"].append(
                {
                    "path": str(mapping),
                    "name": mapping.name,
                    "role": "representative_mapping",
                }
            )

        self._json(result)

    def _api_release(self, query: dict[str, list[str]]) -> None:
        resolved = self.state.resolve(
            config_path=_single(query, "config"),
            snapshot=_single(query, "snapshot"),
        )

        plan = plan_pipeline(resolved, repo_root=self.state.repo_root)

        release = gold_release_summary(plan)

        payload: dict[str, Any] = {
            "published": bool(release),
            "release": release,
            "resolved_config_sha256": resolved.sha256,
            "scientific_config_sha256": resolved.scientific_sha256,
            "snapshot": resolved.get("snapshot.snapshot_id"),
            "snapshot_display": format_snapshot_id(
                str(resolved.get("snapshot.snapshot_id") or "")
            ),
            "brain_threshold_angstrom": resolved.get(
                "brain_filter.threshold_angstrom"
            ),
            "near_duplicate_threshold_angstrom": resolved.get(
                "duplicate_search.near_duplicate_threshold_angstrom"
            ),
            "representative_policy": resolved.get(
                "representative_selection.policy_name"
            ),
        }

        if release:
            explorer_summary: dict[str, Any] = {}

            try:
                explorer_summary = self.state.explorer(resolved).summary()
            except DuplicateQueryError:
                explorer_summary = {}

            payload["pair_counts"] = {
                key: explorer_summary.get(key)
                for key in (
                    "near_duplicate_pairs",
                    "exact_duplicate_pairs",
                    "nonzero_near_duplicate_pairs",
                    "non_near_duplicate_pairs",
                    "total_tested_pairs",
                )
            }

        self._json(payload)

    def _api_runs(self, query: dict[str, list[str]]) -> None:
        resolved = self.state.resolve(config_path=_single(query, "config"))
        paths_run_root = resolved.get("storage.run_root") or "outputs/runs"

        run_root = Path(paths_run_root)

        if not run_root.is_absolute():
            run_root = self.state.repo_root / run_root

        self._json({"run_root": str(run_root), "runs": list_runs(run_root)})

    def _api_run_detail(self, run_id: str) -> None:
        resolved = self.state.resolve()
        run_root = Path(resolved.get("storage.run_root") or "outputs/runs")

        if not run_root.is_absolute():
            run_root = self.state.repo_root / run_root

        directory = run_root / unquote(run_id)

        if not directory.is_dir():
            self._error("No such run", HTTPStatus.NOT_FOUND)
            return

        self._json(RunProvenance.load(directory).record)

    # -- historical run inspection (strictly read-only) -----------------

    def _run_directory(self, run_id: str) -> Path:
        resolved = self.state.resolve()
        run_root = Path(resolved.get("storage.run_root") or "outputs/runs")

        if not run_root.is_absolute():
            run_root = self.state.repo_root / run_root

        directory = run_root / unquote(run_id)

        if not directory.is_dir():
            raise UIError(f"No such run: {run_id}")

        return directory

    def _load_record(self, run_id: str) -> dict[str, Any]:
        """Read one run's record. Never writes, never re-resolves anything."""

        record_path = self._run_directory(run_id) / "run.json"

        if not record_path.is_file():
            raise UIError(f"Run has no run.json: {run_id}")

        return json.loads(record_path.read_text(encoding="utf-8"))

    def _api_run_timeline(self, run_id: str) -> None:
        record = self._load_record(run_id)

        self._json(
            {
                "run_id": record.get("run_id", run_id),
                "created_at": record.get("created_at"),
                "status": record.get("status"),
                "snapshot": record.get("snapshot") or {},
                "resolved_config_sha256": record.get("resolved_config_sha256"),
                "scientific_config_sha256": record.get(
                    "scientific_config_sha256"
                ),
                "git": record.get("git") or {},
                "environment": record.get("environment") or {},
                "runtime": record.get("runtime") or {},
                "config_overrides": record.get("config_overrides") or [],
                "config_file": record.get("config_file"),
                "run_directory": record.get("run_directory"),
                "release": record.get("release") or {},
                # Always canonical scientific order, never alphabetical.
                "timeline": run_timeline(record),
            }
        )

    def _api_run_stage(
        self,
        run_id: str,
        canonical_key: str,
        query: dict[str, list[str]],
    ) -> None:
        record = self._load_record(run_id)

        try:
            detail = stage_detail(
                record,
                unquote(canonical_key),
                repo_root=self.state.repo_root,
                list_artefacts=not _flag(query, "no_artefacts"),
            )
        except KeyError as exc:
            self._error(str(exc), HTTPStatus.NOT_FOUND)
            return

        detail["run_id"] = record.get("run_id", run_id)
        detail["duplicate_navigation"] = duplicate_navigation(
            unquote(canonical_key)
        )
        detail["read_only"] = True

        self._json(detail)

    # -- artefact access (allowlisted, read-only) -----------------------

    def _scene_cache_root(self) -> Path:
        return self.state.scene_cache_root(self.state.resolve())

    def _allowed_artefact_roots(self) -> list[Path]:
        """The only directories the server will ever read a file from.

        Deliberately NOT the repository root: this is a provenance viewer, not
        a filesystem browser. Only pipeline outputs, releases, run provenance,
        the durable/hot snapshot stores and the Mol* report assets are
        reachable.
        """

        resolved = self.state.resolve()
        repo = self.state.repo_root.resolve()

        roots: list[Path] = []

        for dotted in (
            "storage.output_root",
            "storage.release_root",
            "storage.run_root",
            "storage.durable_snapshot_root",
            "storage.hot_cache_root",
        ):
            value = resolved.get(dotted)

            if not value:
                continue

            candidate = Path(value)

            if not candidate.is_absolute():
                candidate = repo / candidate

            roots.append(candidate)

        # Prepared Mol* assets and generated scenes.
        roots.append(repo / "reports" / "molstar_exact_duplicate_examples")
        roots.append(self._scene_cache_root())

        allowed: list[Path] = []

        for root in roots:
            try:
                allowed.append(root.resolve())
            except OSError:  # pragma: no cover - unreadable mount
                continue

        return allowed

    def _resolve_artefact_path(self, raw: str) -> Path:
        """Resolve a requested path, or refuse it.

        ``Path.resolve()`` collapses ``..`` and follows symlinks before the
        containment check, so neither traversal nor a symlink pointing outside
        an allowed root can escape.
        """

        if not raw:
            raise UIError("An artefact path is required")

        target = Path(raw)

        if not target.is_absolute():
            target = self.state.repo_root / target

        try:
            target = target.resolve()
        except OSError as exc:
            raise UIError(f"Could not resolve artefact path: {exc}") from exc

        allowed = self._allowed_artefact_roots()

        if not any(
            target == root or root in target.parents for root in allowed
        ):
            raise ArtefactForbidden(
                "Artefact is outside the directories this viewer may read."
            )

        return target

    def _api_artefact(self, query: dict[str, list[str]]) -> None:
        """Metadata plus a bounded preview of one artefact."""

        try:
            target = self._resolve_artefact_path(_single(query, "path") or "")
        except ArtefactForbidden as exc:
            self._error(str(exc), HTTPStatus.FORBIDDEN)
            return

        if not target.is_file():
            self._error(f"No such artefact: {target}", HTTPStatus.NOT_FOUND)
            return

        payload = artefact_module.describe(target)
        payload["preview_kind"] = artefact_module.preview_kind(target)
        payload["download_url"] = (
            "/api/artefact/download?path=" + quote(str(target))
        )
        payload["provenance"] = self._artefact_provenance(target)

        if payload["preview_kind"] == "parquet":
            try:
                payload["schema"] = artefact_module.parquet_schema(target)
            except artefact_module.ArtefactError as exc:
                self._error(str(exc))
                return
        else:
            try:
                payload["preview"] = artefact_module.preview(
                    target, limit=_int(query, "limit") or 50
                )["preview"]
            except artefact_module.ArtefactError as exc:
                self._error(str(exc))
                return

        self._json(payload)

    def _api_artefact_table(self, query: dict[str, list[str]]) -> None:
        """One bounded page of a tabular artefact."""

        try:
            target = self._resolve_artefact_path(_single(query, "path") or "")
        except ArtefactForbidden as exc:
            self._error(str(exc), HTTPStatus.FORBIDDEN)
            return

        if not target.is_file():
            self._error(f"No such artefact: {target}", HTTPStatus.NOT_FOUND)
            return

        descending = _flag(query, "descending")

        try:
            payload = artefact_module.table_page(
                target,
                page=_int(query, "page") or 1,
                page_size=_int(query, "page_size")
                or artefact_module.DEFAULT_PAGE_SIZE,
                search=_single(query, "search"),
                sort_by=_single(query, "sort_by"),
                descending=descending,
            )
        except artefact_module.ArtefactError as exc:
            self._error(str(exc))
            return

        payload["path"] = str(target)
        payload["name"] = target.name

        self._json(payload)

    def _api_artefact_download(self, query: dict[str, list[str]]) -> None:
        """Serve the original artefact bytes, unmodified."""

        try:
            target = self._resolve_artefact_path(_single(query, "path") or "")
        except ArtefactForbidden as exc:
            self._error(str(exc), HTTPStatus.FORBIDDEN)
            return

        if not target.is_file():
            self._error(f"No such artefact: {target}", HTTPStatus.NOT_FOUND)
            return

        size = target.stat().st_size

        if size > MAX_DOWNLOAD_BYTES:
            self._error(
                f"Artefact is {size} bytes, above the "
                f"{MAX_DOWNLOAD_BYTES}-byte download limit. Copy it from "
                f"{target} directly instead.",
                HTTPStatus.REQUEST_ENTITY_TOO_LARGE,
            )
            return

        # The bytes on disk, byte for byte. Never a re-encoded export.
        body = target.read_bytes()

        content_type = (
            mimetypes.guess_type(target.name)[0] or "application/octet-stream"
        )

        self._send(
            HTTPStatus.OK,
            body,
            content_type,
            extra_headers={
                "Content-Disposition": (
                    f'attachment; filename="{target.name}"'
                ),
                "X-Artefact-Path": str(target),
            },
        )

    def _artefact_provenance(self, target: Path) -> dict[str, Any]:
        """Best-effort provenance for an artefact, from records only.

        Reads nothing but existing summaries and run records, and never
        writes. Absent facts are reported as null rather than invented.
        """

        provenance: dict[str, Any] = {
            "run_id": None,
            "stage": None,
            "snapshot_id": None,
            "resolved_config_sha256": None,
            "scientific_config_sha256": None,
            "producer": None,
            "validation": None,
        }

        # A stage summary sitting beside the artefact names the snapshot and
        # the protocol the artefact belongs to.
        for parent in list(target.parents)[:3]:
            summary = parent / "global_summary.json"

            if summary.is_file():
                try:
                    loaded = json.loads(summary.read_text(encoding="utf-8"))
                except (OSError, ValueError):
                    continue

                provenance["snapshot_id"] = loaded.get("snapshot")
                provenance["stage"] = loaded.get("summary_schema_name")
                break

        # A release manifest names the release identity.
        for parent in list(target.parents)[:3]:
            manifest = parent / "release_manifest.json"

            if manifest.is_file():
                try:
                    loaded = json.loads(manifest.read_text(encoding="utf-8"))
                except (OSError, ValueError):
                    continue

                provenance["snapshot_id"] = loaded.get("snapshot")
                provenance["stage"] = "Stage 14c — Final Gold release"
                provenance["release_name"] = loaded.get("release_name")
                break

        # A run directory names the run and both configuration hashes.
        for parent in list(target.parents)[:4]:
            record_path = parent / "run.json"

            if record_path.is_file():
                try:
                    record = json.loads(
                        record_path.read_text(encoding="utf-8")
                    )
                except (OSError, ValueError):
                    continue

                provenance["run_id"] = record.get("run_id")
                provenance["resolved_config_sha256"] = record.get(
                    "resolved_config_sha256"
                )
                provenance["scientific_config_sha256"] = record.get(
                    "scientific_config_sha256"
                )
                provenance["snapshot_id"] = (
                    record.get("snapshot") or {}
                ).get("snapshot_id") or provenance["snapshot_id"]
                break

        return provenance

    # -- Mol* on demand -------------------------------------------------

    def _pair_from_query(self, query: dict[str, list[str]]) -> PairRequest:
        def _get(name: str) -> str:
            value = _single(query, name)

            if not value:
                raise UIError(f"Missing pair parameter: {name}")

            return value

        return PairRequest(
            pdb_id_a=_get("pdb_id_a"),
            chain_a=_get("chain_a"),
            pdb_id_b=_get("pdb_id_b"),
            chain_b=_get("chain_b"),
            snapshot_id=_single(query, "snapshot"),
            run_id=_single(query, "run_id"),
            model_id=_int(query, "model_id") or 1,
            chain_length=_int(query, "chain_length"),
            d_bri_units=_int(query, "d_bri"),
            classification=_single(query, "classification"),
            relationship=_single(query, "relationship"),
            representative=_single(query, "representative"),
        )

    def _api_pair_availability(self, query: dict[str, list[str]]) -> None:
        """Can this pair be visualised, and if not, precisely why not."""

        resolved = self.state.resolve(snapshot=_single(query, "snapshot"))
        request = self._pair_from_query(query)

        service = self.state.molstar(
            resolved, snapshot_id=request.snapshot_id
        )

        payload = service.availability(request)
        payload["pair"] = request.to_dict()

        self._json(payload)

    def _api_pair_scene(self, query: dict[str, list[str]]) -> None:
        """Generate (or reuse) a Mol* scene for one recorded pair."""

        resolved = self.state.resolve(snapshot=_single(query, "snapshot"))
        request = self._pair_from_query(query)
        view = _single(query, "view") or "side_by_side"

        service = self.state.molstar(
            resolved, snapshot_id=request.snapshot_id
        )

        try:
            payload = service.scene(request, view)
        except MolstarServiceError as exc:
            self._error(str(exc), HTTPStatus.NOT_FOUND)
            return

        # Structures are served from the allowlisted roots by /structures/.
        payload["structure_urls"] = {
            source["pdb_id"]: (
                "/api/structure?path=" + quote(source["path"])
            )
            for source in (
                payload["provenance"]["structure_a"],
                payload["provenance"]["structure_b"],
            )
        }

        self._json(payload)

    def _api_structure(self, query: dict[str, list[str]]) -> None:
        """Serve one allowlisted structure file to the viewer."""

        try:
            target = self._resolve_artefact_path(_single(query, "path") or "")
        except ArtefactForbidden as exc:
            self._error(str(exc), HTTPStatus.FORBIDDEN)
            return

        if not target.is_file() or target.suffix.lower() not in {
            ".cif",
            ".bcif",
            ".pdb",
            ".mmcif",
        }:
            self._error("Not a servable structure", HTTPStatus.NOT_FOUND)
            return

        self._send(
            HTTPStatus.OK,
            target.read_bytes(),
            "chemical/x-mmcif",
            extra_headers={"X-Structure-Path": str(target)},
        )

    def _api_scenes(self) -> None:
        self._json({"scenes": _scene_index(self.state.repo_root)})

    # -- POST endpoints -------------------------------------------------

    def _api_resolve(self, body: dict[str, Any]) -> None:
        resolved = self.state.resolve(
            config_path=body.get("config_path"),
            overrides=body.get("overrides"),
            snapshot=body.get("snapshot"),
        )

        payload: dict[str, Any] = {
            "resolved": resolved.to_dict(),
            "resolved_config_sha256": resolved.sha256,
            "scientific_config_sha256": resolved.scientific_sha256,
            "resolved_config_yaml": resolved.to_yaml(),
            "sources": resolved.sources,
            "layers": list(resolved.layers),
            "scientific_summary": resolved.scientific_summary(),
            "near_duplicate_threshold_mA": (
                resolved.near_duplicate_threshold_mA
            ),
            "brain_threshold_mA": resolved.brain_threshold_mA,
            "representation_precision_angstrom": resolved.get(
                "bri.representation_precision_angstrom"
            ),
            "representation_unit": representation_unit_label(
                float(resolved.get("bri.representation_precision_angstrom"))
            ),
            "precision_is_implemented": precision_is_implemented(
                resolved.data
            ),
            "implemented_precision_angstrom": (
                IMPLEMENTED_REPRESENTATION_PRECISION_ANGSTROM
            ),
        }

        if resolved.get("snapshot.snapshot_id"):
            paths = self.state.paths(resolved)

            layout = SnapshotStoreLayout.from_config(
                resolved, repo_root=self.state.repo_root
            )

            payload["snapshot_status"] = snapshot_status(
                str(resolved.get("snapshot.snapshot_id")),
                layout,
                remote_available=True,
            )

            payload["paths"] = {
                "output_root": str(paths.output_root),
                "release_root": str(paths.release_root),
                "run_root": str(paths.run_root),
                "release": paths.release,
            }

        self._json(payload)

    def _api_plan(self, body: dict[str, Any]) -> None:
        resolved = self.state.resolve(
            config_path=body.get("config_path"),
            overrides=body.get("overrides"),
            snapshot=body.get("snapshot"),
        )

        if not resolved.get("snapshot.snapshot_id"):
            self._error(
                "Select a snapshot before planning; a run must be pinned to a "
                "concrete snapshot identity."
            )
            return

        plan = plan_pipeline(resolved, repo_root=self.state.repo_root)

        self._json(plan.to_dict())

    def _api_run(self, body: dict[str, Any]) -> None:
        """Create a run: freeze identity and write provenance.

        The UI deliberately creates the run and records provenance, then
        reports the exact commands.  Submitting long Slurm work is done from
        the CLI on a login node, so a browser tab can never leave a partially
        submitted array behind.
        """

        resolved = self.state.resolve(
            config_path=body.get("config_path"),
            overrides=body.get("overrides"),
            snapshot=body.get("snapshot"),
        )

        if not resolved.get("snapshot.snapshot_id"):
            self._error("Select a snapshot before starting a run.")
            return

        paths = self.state.paths(resolved)
        plan = plan_pipeline(resolved, repo_root=self.state.repo_root)

        provenance = RunProvenance.create(
            resolved=resolved,
            run_root=paths.run_root,
            repo_root=self.state.repo_root,
            snapshot={
                "snapshot_id": resolved.get("snapshot.snapshot_id"),
                "display": format_snapshot_id(
                    str(resolved.get("snapshot.snapshot_id"))
                ),
                "selection_mode": resolved.get(
                    "snapshot.resolved_selection_mode"
                ),
            },
            invocation={"origin": "ui", "executor": "deferred"},
        )

        record_plan_in_provenance(plan, provenance)
        provenance.set_status("planned")
        provenance.flush()

        from pdbclean.cli import stage_command

        commands = []

        for observation in plan.to_run:
            try:
                argv = stage_command(
                    observation.stage.stage_id,
                    resolved,
                    paths,
                )
            except PipelineError as exc:
                argv = None
                commands.append(
                    {
                        "stage_id": observation.stage.stage_id,
                        "error": str(exc),
                    }
                )
                continue

            if argv:
                commands.append(
                    {
                        "stage_id": observation.stage.stage_id,
                        "argv": argv,
                    }
                )

        self._json(
            {
                "run_id": provenance.run_id,
                "run_directory": str(provenance.run_dir),
                "resolved_config_sha256": resolved.sha256,
                "scientific_config_sha256": resolved.scientific_sha256,
                "plan": plan.to_dict(),
                "commands": commands,
                "cli_equivalent": (
                    "pdbclean run"
                    + (
                        f" --config {body['config_path']}"
                        if body.get("config_path")
                        else ""
                    )
                    + f" --snapshot {resolved.get('snapshot.snapshot_id')}"
                    + " --executor slurm"
                ),
            }
        )

    # -- static and structure files -------------------------------------

    def _serve_static(self, route: str) -> None:
        if route in ("", "/"):
            route = "/index.html"

        relative = posixpath.normpath(unquote(route)).lstrip("/")

        if relative.startswith(".."):
            self._error("Invalid path", HTTPStatus.FORBIDDEN)
            return

        candidate = (STATIC_DIR / relative).resolve()

        try:
            candidate.relative_to(STATIC_DIR.resolve())
        except ValueError:
            self._error("Invalid path", HTTPStatus.FORBIDDEN)
            return

        if not candidate.is_file():
            self._error("Not found", HTTPStatus.NOT_FOUND)
            return

        content_type, _ = mimetypes.guess_type(str(candidate))

        self._send(
            HTTPStatus.OK,
            candidate.read_bytes(),
            content_type or "application/octet-stream",
        )

    def _serve_structure(self, route: str) -> None:
        relative = posixpath.normpath(
            unquote(route[len("/structures/"):])
        ).lstrip("/")

        if relative.startswith(".."):
            self._error("Invalid path", HTTPStatus.FORBIDDEN)
            return

        repo_root = self.state.repo_root.resolve()

        for root in STRUCTURE_ROOTS:
            base = (repo_root / root).resolve()
            candidate = (base / relative).resolve()

            try:
                candidate.relative_to(base)
            except ValueError:
                continue

            if not candidate.is_file():
                continue

            if candidate.suffix.lower() not in STRUCTURE_SUFFIXES:
                self._error("Unsupported file type", HTTPStatus.FORBIDDEN)
                return

            content_type, _ = mimetypes.guess_type(str(candidate))

            self._send(
                HTTPStatus.OK,
                candidate.read_bytes(),
                content_type or "chemical/x-cif",
            )
            return

        self._error("Not found", HTTPStatus.NOT_FOUND)


def _single(query: dict[str, list[str]], key: str) -> str | None:
    values = query.get(key)

    if not values:
        return None

    value = values[0].strip()

    return value or None


def _flag(query: dict[str, list[str]], key: str) -> bool:
    value = _single(query, key)

    return value is not None and value.lower() in {"1", "true", "yes", "on"}


def _int(query: dict[str, list[str]], key: str) -> int | None:
    value = _single(query, key)

    if value is None:
        return None

    try:
        return int(value)
    except ValueError as exc:
        raise UIError(f"{key} must be an integer") from exc


def _scene_index(repo_root: Path) -> list[dict[str, Any]]:
    """Index the prepared Mol* scenes already present in the repository.

    The existing ``reports/molstar_exact_duplicate_examples`` work is reused as
    is.  Nothing is regenerated and nothing is moved.
    """

    base = repo_root / "reports" / "molstar_exact_duplicate_examples"

    metrics_path = base / "metrics.json"

    if not metrics_path.is_file():
        return []

    try:
        metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return []

    scenes: list[dict[str, Any]] = []

    for pair_key, payload in sorted(metrics.items()):
        reference = str(payload.get("reference", ""))
        moving = str(payload.get("moving", ""))

        if ":" not in reference or ":" not in moving:
            continue

        pdb_a, chain_a = reference.split(":", 1)
        pdb_b, chain_b = moving.split(":", 1)

        views = []

        for path in sorted(base.glob(f"{pair_key}_*.mvsj")):
            label = path.stem[len(pair_key) + 1:].replace("_", " ")

            views.append(
                {
                    "label": label,
                    "url": f"/structures/{path.name}",
                }
            )

        if not views:
            continue

        scenes.append(
            {
                "key": pair_key,
                "pdb_id_a": pdb_a.lower(),
                "chain_a": chain_a,
                "pdb_id_b": pdb_b.lower(),
                "chain_b": chain_b,
                "d_bri_mA": payload.get("stage10_d_bri_mA"),
                "d_bri_angstrom": payload.get("stage10_d_bri_A"),
                "chain_length": payload.get("retained_residue_count"),
                "backbone_atom_count": payload.get("backbone_atom_count"),
                "is_zero_duplicate": payload.get(
                    "stage10_is_zero_duplicate"
                ),
                "views": views,
            }
        )

    return scenes


def serve(
    *,
    host: str = "127.0.0.1",
    port: int = 8765,
    repo_root: Path,
    config_path: str | None = None,
    overrides: list[str] | None = None,
    open_browser: bool = True,
) -> int:
    """Run the UI server until interrupted."""

    if not STATIC_DIR.is_dir():
        raise UIError(f"UI assets are missing: {STATIC_DIR}")

    Handler.state = UIState(
        repo_root=Path(repo_root).resolve(),
        config_path=config_path,
        overrides=list(overrides or []),
    )

    server = ThreadingHTTPServer((host, port), Handler)

    url = f"http://{host}:{port}/"

    print("PDBClean UI")
    print(f"  repository : {Handler.state.repo_root}")
    print(f"  config     : {config_path or '(built-in validated defaults)'}")
    print(f"  address    : {url}")
    print("  press Ctrl-C to stop")

    if open_browser:
        try:
            webbrowser.open(url)
        except Exception:
            pass

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nstopping")
    finally:
        server.server_close()

    return 0
