/* PDBClean UI.
 *
 * Presentation only. Every value shown here is produced by the Python backend
 * (pdbclean.runconfig, pdbclean.pipeline, pdbclean.duplicates,
 * pdbclean.run_provenance). No scientific decision is made in this file.
 */

"use strict";

const state = {
    bootstrap: null,
    resolved: null,
    plan: null,
    snapshot: "",
    dupOffset: 0,
    dupMatched: 0,
    scenes: [],
};

/* ------------------------------------------------------------ utilities */

function $(id) {
    return document.getElementById(id);
}

function el(tag, attrs, children) {
    const node = document.createElement(tag);

    for (const key in (attrs || {})) {
        if (key === "class") {
            node.className = attrs[key];
        } else if (key === "text") {
            node.textContent = attrs[key];
        } else if (attrs[key] !== null && attrs[key] !== undefined) {
            node.setAttribute(key, attrs[key]);
        }
    }

    (children || []).forEach(function (child) {
        node.appendChild(
            typeof child === "string"
                ? document.createTextNode(child)
                : child
        );
    });

    return node;
}

function num(value) {
    if (value === null || value === undefined) {
        return "—";
    }

    return Number(value).toLocaleString("en-GB");
}

function toast(message, isError) {
    const box = $("toast");

    box.textContent = message;
    box.style.borderLeftColor = isError
        ? "var(--fail)"
        : "var(--accent)";
    box.hidden = false;

    clearTimeout(toast._timer);
    toast._timer = setTimeout(function () {
        box.hidden = true;
    }, 6000);
}

async function api(path, options) {
    const response = await fetch(path, options);
    const payload = await response.json();

    if (!response.ok || payload.error) {
        throw new Error(payload.error || response.statusText);
    }

    return payload;
}

function stateCell(value) {
    return el("span", {
        class: "state state-" + (value || "pending"),
        text: (value || "pending").replace(/_/g, " "),
    });
}

/* ---------------------------------------------------------------- theme
 *
 * Purely a viewing preference. The theme is stored in localStorage and
 * applied by setting data-theme on <html>. It is never sent to the backend,
 * never enters the resolved run configuration, and cannot affect
 * resolved_config_sha256 or scientific_config_sha256.
 */

const THEME_KEY = "pdbclean.theme";

function storedTheme() {
    try {
        return window.localStorage.getItem(THEME_KEY);
    } catch (error) {
        return null;   /* private mode: fall back to the OS preference */
    }
}

function systemPrefersDark() {
    return (
        window.matchMedia
        && window.matchMedia("(prefers-color-scheme: dark)").matches
    );
}

function activeTheme() {
    return storedTheme() || (systemPrefersDark() ? "dark" : "light");
}

function applyTheme(theme) {
    if (theme) {
        document.documentElement.setAttribute("data-theme", theme);
    } else {
        document.documentElement.removeAttribute("data-theme");
    }

    const button = $("theme-toggle");

    if (button) {
        const next = theme === "dark" ? "light" : "dark";

        button.textContent = theme === "dark" ? "light theme" : "dark theme";
        button.setAttribute("aria-label", "Switch to the " + next + " theme");
    }
}

function initTheme() {
    /* Applied before first paint by the inline bootstrap in index.html; this
     * only syncs the button label and wires the toggle. */
    applyTheme(activeTheme());

    const button = $("theme-toggle");

    if (!button) {
        return;
    }

    button.addEventListener("click", function () {
        const next = activeTheme() === "dark" ? "light" : "dark";

        try {
            window.localStorage.setItem(THEME_KEY, next);
        } catch (error) {
            /* Not persisting is acceptable; the switch still applies. */
        }

        /* Only swaps CSS custom properties. No state is reloaded or lost. */
        applyTheme(next);
    });
}

/* --------------------------------------------------------------- routing */

function showView(name) {
    document.querySelectorAll("main > section").forEach(function (node) {
        node.hidden = node.id !== "view-" + name;
    });

    document.querySelectorAll("nav button").forEach(function (button) {
        button.setAttribute(
            "aria-current",
            button.dataset.view === name ? "true" : "false"
        );
    });

    if (name === "pipeline" && !state.plan) {
        loadPlan();
    }

    if (name === "duplicates" && !$("dup-table").tBodies[0].rows.length) {
        searchDuplicates(0);
    }

    if (name === "release") {
        loadRelease();
    }

    if (name === "runs") {
        loadRuns();
    }
}

/* --------------------------------------------------------- configuration */

function currentOverrides() {
    const overrides = {};

    function put(key, raw, cast) {
        const value = (raw || "").trim();

        if (value === "") {
            return;
        }

        overrides[key] = cast === "number" ? Number(value) : value;
    }

    put("selection.models.model_id", $("cfg-model").value, "number");
    put(
        "quality_rules.backbone_distance.minimum_distance_angstrom",
        $("cfg-q005").value,
        "number"
    );
    put(
        "post_cleaning_geometric_validation.minimum_triangle_angle_degrees",
        $("cfg-angle").value,
        "number"
    );
    put(
        "bri.representation_precision_angstrom",
        $("cfg-precision").value,
        "number"
    );
    put("brain_filter.threshold_angstrom", $("cfg-brain").value, "number");
    put(
        "duplicate_search.near_duplicate_threshold_angstrom",
        $("cfg-near").value,
        "number"
    );

    return overrides;
}

function selectedSnapshot() {
    const manual = $("cfg-snapshot-manual").value.trim();

    if (manual) {
        return manual;
    }

    return $("cfg-snapshot").value || "";
}

function requestBody() {
    return {
        config_path: $("cfg-profile").value || null,
        overrides: currentOverrides(),
        snapshot: selectedSnapshot() || null,
    };
}

async function resolveConfiguration(quiet) {
    $("config-error").innerHTML = "";

    try {
        const payload = await api("/api/config/resolve", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify(requestBody()),
        });

        state.resolved = payload;
        state.plan = null;
        state.snapshot = payload.resolved.snapshot.snapshot_id || "";

        renderResolved(payload);

        if (!quiet) {
            toast("Configuration resolved.");
        }

        return payload;
    } catch (error) {
        $("config-error").appendChild(
            el("div", { class: "fail-block", text: String(error.message) })
        );

        if (!quiet) {
            toast(error.message, true);
        }

        return null;
    }
}

function renderResolved(payload) {
    const body = $("resolved-table").tBodies[0];

    body.innerHTML = "";

    payload.scientific_summary.forEach(function (row) {
        const origin = row.source || "builtin_default";

        body.appendChild(
            el("tr", {}, [
                el("td", { text: row.label }),
                el("td", { class: "mono", text: String(row.value) }),
                el("td", {
                    class: origin.indexOf("builtin") === 0 ? "note" : "origin",
                    text: origin,
                }),
            ])
        );
    });

    const identity = $("resolved-identity");

    identity.innerHTML = "";

    const facts = [
        ["Resolved config SHA256", payload.resolved_config_sha256],
        ["Scientific SHA256", payload.scientific_config_sha256],
        ["Brain threshold (mÅ)", String(payload.brain_threshold_mA)],
        [
            "Complete-BRI threshold (mÅ)",
            String(payload.near_duplicate_threshold_mA),
        ],
        [
            "BRI representation precision",
            payload.representation_precision_angstrom + " Å  (1 unit = 1 "
                + payload.representation_unit + ")",
        ],
        ["Configuration layers", payload.layers.join("  →  ")],
    ];

    if (payload.snapshot_status) {
        const snap = payload.snapshot_status;

        facts.push(["Snapshot availability", snap.availability]);
        facts.push([
            "Reproducible without cache",
            snap.reproducible_without_cache ? "yes" : "not yet preserved",
        ]);
    }

    if (payload.paths) {
        facts.push(["Release name", payload.paths.release]);
        facts.push(["Output root", payload.paths.output_root]);
        facts.push(["Run root", payload.paths.run_root]);
    }

    facts.forEach(function (pair) {
        identity.appendChild(el("dt", { text: pair[0] }));
        identity.appendChild(el("dd", { text: pair[1] }));
    });

    renderPrecisionNotice(payload);

    $("resolved-yaml").textContent = payload.resolved_config_yaml;

    $("fact-snapshot").textContent =
        payload.resolved.snapshot.snapshot_id || "not selected";
    $("fact-threshold").textContent =
        payload.near_duplicate_threshold_mA + " mÅ";
    $("fact-sha").textContent =
        payload.scientific_config_sha256.slice(0, 12);
    $("fact-defaults").textContent = payload.resolved.defaults_version;
}

/* Make an experimental precision unmistakable, without implying it is
 * equally validated. */
function renderPrecisionNotice(payload) {
    const holder = $("precision-note");

    if (!holder) {
        return;
    }

    holder.innerHTML = "";

    const precision =
        payload.resolved.bri.representation_precision_angstrom;
    const validated =
        (state.bootstrap && state.bootstrap.defaults.bri
            .representation_precision_angstrom) || 0.001;

    if (Number(precision) === Number(validated)) {
        return;
    }

    holder.appendChild(
        el("div", {
            class: "fail-block",
            text:
                "Experimental configuration: BRI representation precision is "
                + precision + " A, not the validated " + validated + " A. "
                + "This is a distinct scientific configuration with its own "
                + "scientific identity; it cannot reuse the validated "
                + "release, and the production stages do not implement this "
                + "grid.",
        })
    );
}

function fillFormFromDefaults(defaults) {
    $("cfg-model").value = defaults.selection.models.model_id;
    $("cfg-q005").value =
        defaults.quality_rules.backbone_distance.minimum_distance_angstrom;
    $("cfg-angle").value =
        defaults.post_cleaning_geometric_validation
            .minimum_triangle_angle_degrees;
    $("cfg-precision").value = defaults.bri.representation_precision_angstrom;
    $("cfg-brain").value = defaults.brain_filter.threshold_angstrom;
    $("cfg-near").value =
        defaults.duplicate_search.near_duplicate_threshold_angstrom;
}

async function listSnapshots() {
    $("snapshot-note").textContent = "querying the archive…";

    try {
        const payload = await api(
            "/api/snapshots?limit=40"
            + ($("cfg-profile").value
                ? "&config=" + encodeURIComponent($("cfg-profile").value)
                : "")
        );

        const select = $("cfg-snapshot");
        const chosen = select.value;

        select.innerHTML = "";
        select.appendChild(
            el("option", {
                value: "",
                text: "latest complete snapshot (default)",
            })
        );

        payload.snapshots.forEach(function (entry) {
            select.appendChild(
                el("option", {
                    value: entry.snapshot_id,
                    text:
                        entry.display
                        + (entry.is_latest ? "  (latest discovered)" : ""),
                })
            );
        });

        select.value = chosen;

        $("snapshot-note").textContent =
            payload.error
                ? payload.error
                : payload.snapshots.length + " snapshots discovered in "
                  + payload.bucket_url;
    } catch (error) {
        $("snapshot-note").textContent = error.message;
    }
}

/* ---------------------------------------------------------------- plan */

async function loadPlan() {
    const status = $("plan-status");

    status.innerHTML = "";
    $("plan-table").tBodies[0].innerHTML = "";

    try {
        const plan = await api("/api/plan", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify(requestBody()),
        });

        state.plan = plan;
        renderPlan(plan);
    } catch (error) {
        status.appendChild(
            el("div", { class: "fail-block", text: error.message })
        );
    }
}

function renderLifecycle(plan) {
    const holder = $("lifecycle");

    holder.innerHTML = "";

    const layers = state.bootstrap.layers;

    layers.forEach(function (layer) {
        const stages = plan.stages.filter(function (stage) {
            return stage.layer === layer.id;
        });

        const passed = stages.filter(function (stage) {
            return (
                stage.validation === "validation_pass"
                || stage.validation === "not_applicable"
            );
        }).length;

        holder.appendChild(
            el("div", {}, [
                el("div", { class: "name", text: layer.label }),
                el("div", {
                    class: "detail",
                    text: passed + " / " + stages.length + " gates passed",
                }),
            ])
        );
    });
}

function renderPlan(plan) {
    renderLifecycle(plan);

    const status = $("plan-status");

    status.innerHTML = "";

    const blocked = plan.stages.filter(function (stage) {
        return stage.action === "blocked";
    });

    const failing = plan.stages.filter(function (stage) {
        return stage.incompatibilities.length > 0;
    });

    if (failing.length) {
        const box = el("div", { class: "warn-block" }, [
            el("div", {
                text:
                    failing.length
                    + " stage(s) have outputs produced under a different "
                    + "scientific configuration. They will not be reused.",
            }),
        ]);

        failing.forEach(function (stage) {
            stage.incompatibilities.forEach(function (issue) {
                box.appendChild(
                    el("div", {
                        class: "note",
                        text:
                            stage.stage_id + ": " + issue.summary_key
                            + " is " + JSON.stringify(issue.observed)
                            + ", configuration requires "
                            + JSON.stringify(issue.expected),
                    })
                );
            });
        });

        status.appendChild(box);
    }

    if (blocked.length) {
        status.appendChild(
            el("div", {
                class: "warn-block",
                text:
                    blocked.length
                    + " downstream stage(s) are blocked until upstream "
                    + "validation passes.",
            })
        );
    }

    const body = $("plan-table").tBodies[0];

    body.innerHTML = "";

    plan.stages.forEach(function (stage) {
        const row = el("tr", { class: "stage-row" }, [
            el("td", { class: "canonical" }, [
                el("span", { class: "canonical-label", text: stage.canonical_stage }),
            ]),
            el("td", {}, [
                el("button", {
                    class: "link",
                    text: stage.title,
                    "data-stage": stage.stage_id,
                }),
            ]),
            el("td", {}, [
                el("span", { class: "layer-tag", text: stage.layer }),
            ]),
            el("td", {}, [stateCell(stage.status)]),
            el("td", {}, [stateCell(stage.validation)]),
            el("td", { text: stage.action.replace(/_/g, " ") }),
            el("td", { class: "num", text: num(stage.input_count) }),
            el("td", { class: "num", text: num(stage.output_count) }),
            el("td", {
                class: "num",
                text: stage.primary_output_bytes
                    ? num(stage.primary_output_bytes)
                    : "—",
            }),
            el("td", {
                class: "mono",
                text: stage.slurm_job_ids.join(", ") || "—",
            }),
        ]);

        body.appendChild(row);

        row.querySelector("button.link").addEventListener(
            "click",
            function () {
                toggleStageDetail(row, stage);
            }
        );
    });
}

function toggleStageDetail(row, stage) {
    const next = row.nextElementSibling;

    if (next && next.classList.contains("detail-row")) {
        next.remove();
        return;
    }

    const cell = el("td", { colspan: "10" });

    cell.appendChild(el("div", { text: stage.purpose }));

    if (stage.validation_description) {
        cell.appendChild(
            el("div", { text: "Gate: " + stage.validation_description })
        );
    }

    const kv = el("dl", { class: "kv" });

    function add(label, value) {
        if (value === null || value === undefined || value === "") {
            return;
        }

        kv.appendChild(el("dt", { text: label }));
        kv.appendChild(el("dd", { text: String(value) }));
    }

    add("Entry point", stage.entry_point);
    add("Output path", stage.output_path);
    add("Manifest", stage.manifest_path);
    add("Summary", stage.summary_path);
    add("Checksum", stage.primary_output_sha256);
    add("Depends on", stage.depends_on.join(", "));

    Object.keys(stage.scientific_parameters).forEach(function (key) {
        add(key, stage.scientific_parameters[key]);
    });

    cell.appendChild(kv);

    stage.messages.forEach(function (message) {
        cell.appendChild(el("div", { class: "note", text: message }));
    });

    row.after(el("tr", { class: "detail-row" }, [cell]));
}

/* ---------------------------------------------------- duplicate explorer */

function duplicateQuery(offset) {
    const params = new URLSearchParams();

    if ($("cfg-profile").value) {
        params.set("config", $("cfg-profile").value);
    }

    if (state.snapshot) {
        params.set("snapshot", state.snapshot);
    }

    if ($("dup-pdb").value.trim()) {
        params.set("pdb_id", $("dup-pdb").value.trim());
    }

    if ($("dup-chain").value.trim()) {
        params.set("chain", $("dup-chain").value.trim());
    }

    const klass = $("dup-class").value;

    if (klass === "exact") {
        params.set("exact_only", "1");
    } else if (klass === "nonzero") {
        params.set("nonzero_near_only", "1");
    }

    if ($("dup-rel").value) {
        params.set("relationship", $("dup-rel").value);
    }

    [
        ["dup-minlen", "min_length"],
        ["dup-maxlen", "max_length"],
        ["dup-mind", "min_distance"],
        ["dup-maxd", "max_distance"],
    ].forEach(function (pair) {
        const value = $(pair[0]).value.trim();

        if (value !== "") {
            params.set(pair[1], value);
        }
    });

    params.set("offset", String(offset));
    params.set("limit", $("dup-limit").value);

    return params.toString();
}

async function searchDuplicates(offset) {
    const summary = $("dup-summary");

    summary.innerHTML = "";
    $("dup-table").tBodies[0].innerHTML = "";
    $("dup-range").textContent = "searching…";

    try {
        const payload = await api("/api/duplicates?" + duplicateQuery(offset));

        state.dupOffset = payload.offset;
        state.dupMatched = payload.matched;
        state.scenes = payload.scenes || [];

        renderDuplicateSummary(payload);
        renderDuplicateRows(payload);
    } catch (error) {
        $("dup-range").textContent = "—";
        summary.appendChild(
            el("div", { class: "fail-block", text: error.message })
        );
    }
}

function renderDuplicateSummary(payload) {
    const s = payload.summary || {};
    const holder = $("dup-summary");

    const table = el("table", { class: "compact" }, [
        el("thead", {}, [
            el("tr", {}, [
                el("th", { text: "Population" }),
                el("th", { class: "num", text: "Pairs" }),
            ]),
        ]),
    ]);

    const body = el("tbody");

    [
        ["Total tested pairs", s.total_tested_pairs],
        ["Near duplicates (d ≤ threshold)", s.near_duplicate_pairs],
        ["Exact duplicates (d = 0)", s.exact_duplicate_pairs],
        ["Non-zero near duplicates", s.nonzero_near_duplicate_pairs],
        ["Not near duplicates", s.non_near_duplicate_pairs],
    ].forEach(function (pair) {
        if (pair[1] === undefined || pair[1] === null) {
            return;
        }

        body.appendChild(
            el("tr", {}, [
                el("td", { text: pair[0] }),
                el("td", { class: "num", text: num(pair[1]) }),
            ])
        );
    });

    table.appendChild(body);
    holder.appendChild(table);

    if (s.threshold_mA !== undefined && s.threshold_mA !== null) {
        holder.appendChild(
            el("p", {
                class: "note",
                text:
                    "Threshold recorded by the pipeline: d_bri_mA ≤ "
                    + s.threshold_mA
                    + " (" + s.threshold_angstrom + " Å), inclusive, "
                    + "on complete BRI.",
            })
        );
    }

    if (!payload.has_representative_mapping) {
        holder.appendChild(
            el("p", {
                class: "note",
                text:
                    "Stage 14 has not been published for this configuration, "
                    + "so retained/removed relationships are not shown.",
            })
        );
    }
}

function sceneFor(row) {
    return state.scenes.find(function (scene) {
        return (
            (scene.pdb_id_a === row.pdb_id_a
                && scene.chain_a === row.chain_a
                && scene.pdb_id_b === row.pdb_id_b
                && scene.chain_b === row.chain_b)
            || (scene.pdb_id_a === row.pdb_id_b
                && scene.chain_a === row.chain_b
                && scene.pdb_id_b === row.pdb_id_a
                && scene.chain_b === row.chain_a)
        );
    });
}

function renderDuplicateRows(payload) {
    const body = $("dup-table").tBodies[0];

    body.innerHTML = "";

    payload.rows.forEach(function (row) {
        const scene = sceneFor(row);

        const viewCell = el("td");

        if (scene) {
            const link = el("a", {
                href:
                    "/viewer.html?scene=" + encodeURIComponent(scene.key)
                    + "&a=" + encodeURIComponent(row.pdb_id_a + ":" + row.chain_a)
                    + "&b=" + encodeURIComponent(row.pdb_id_b + ":" + row.chain_b)
                    + "&d=" + encodeURIComponent(String(row.d_bri_mA))
                    + "&cls=" + encodeURIComponent(row.classification)
                    + "&rel=" + encodeURIComponent(row.relationship)
                    + "&rep=" + encodeURIComponent(row.representative || "")
                    + "&sha=" + encodeURIComponent(
                        state.resolved
                            ? state.resolved.scientific_config_sha256
                            : ""
                    ),
                target: "_blank",
                rel: "noopener",
                text: "View pair",
            });

            viewCell.appendChild(link);
        } else {
            viewCell.appendChild(
                el("span", { class: "note", text: "no prepared scene" })
            );
        }

        body.appendChild(
            el("tr", {}, [
                el("td", { class: "mono", text: row.pdb_id_a }),
                el("td", { class: "mono", text: row.chain_a }),
                el("td", { class: "mono", text: row.pdb_id_b }),
                el("td", { class: "mono", text: row.chain_b }),
                el("td", { class: "num", text: String(row.model_a) }),
                el("td", { class: "num", text: num(row.chain_length) }),
                el("td", { class: "num", text: String(row.d_bri_mA) }),
                el("td", {
                    class: "num",
                    text: row.d_bri_angstrom.toFixed(4),
                }),
                el("td", {
                    text:
                        row.classification === "exact_duplicate"
                            ? "exact duplicate"
                            : "non-zero near duplicate",
                }),
                el("td", { text: row.relationship }),
                el("td", {
                    class: "mono",
                    text: row.representative || "—",
                }),
                el("td", {
                    class: "num",
                    text:
                        row.direct_d_bri_mA === null
                            ? "—"
                            : String(row.direct_d_bri_mA),
                }),
                viewCell,
            ])
        );
    });

    const first = payload.matched === 0 ? 0 : payload.offset + 1;
    const last = Math.min(payload.offset + payload.rows.length, payload.matched);

    $("dup-range").textContent =
        num(first) + "–" + num(last) + " of " + num(payload.matched);

    $("dup-prev").disabled = payload.offset <= 0;
    $("dup-next").disabled = last >= payload.matched;
}

/* -------------------------------------------------------------- release */

async function loadRelease() {
    const holder = $("release-body");

    holder.innerHTML = "";

    try {
        const params = new URLSearchParams();

        if ($("cfg-profile").value) {
            params.set("config", $("cfg-profile").value);
        }

        if (state.snapshot) {
            params.set("snapshot", state.snapshot);
        }

        const payload = await api("/api/release?" + params.toString());

        if (!payload.published) {
            holder.appendChild(
                el("div", {
                    class: "warn-block",
                    text:
                        "No Gold release has been published for this "
                        + "configuration. Nothing is shown, because a release "
                        + "page must never display values for a run that has "
                        + "not completed.",
                })
            );
            return;
        }

        const release = payload.release;

        const kv = el("dl", { class: "kv" });

        function add(label, value) {
            if (value === null || value === undefined) {
                return;
            }

            kv.appendChild(el("dt", { text: label }));
            kv.appendChild(el("dd", { text: String(value) }));
        }

        add("Snapshot", payload.snapshot_display);
        add("Release", release.release_name);
        add("Release path", release.release_directory);
        add("Retained manifest", release.retained_manifest);
        add("Canonical input chains", num(release.canonical_input_chain_count));
        add("Retained chains", num(release.retained_chain_count));
        add("Removed chains", num(release.removed_chain_count));
        add("m = 1 chains retained", num(release.m1_retained_chain_count));

        const counts = payload.pair_counts || {};

        add("Near-duplicate pairs", num(counts.near_duplicate_pairs));
        add("Exact pairs", num(counts.exact_duplicate_pairs));
        add("Non-zero near pairs", num(counts.nonzero_near_duplicate_pairs));
        add("Non-near tested pairs", num(counts.non_near_duplicate_pairs));
        add("Total tested pairs", num(counts.total_tested_pairs));

        add("Brain threshold (Å)", payload.brain_threshold_angstrom);
        add(
            "Complete-BRI threshold",
            release.near_duplicate_threshold
        );
        add("Near-duplicate relation", release.near_duplicate_relation);
        add("Distance representation", release.distance_representation);
        add("Representative policy", release.representative_policy);
        add("Policy version", release.representative_policy_version);
        add(
            "Direct edge required for removal",
            release.every_removed_chain_has_direct_representative_edge
        );
        add(
            "Components treated as equivalence",
            release.connectedness_treated_as_duplicate_equivalence
        );
        add(
            "Automatic transitive removal",
            release.automatic_transitive_removal
        );
        add("Old-snapshot comparison used", release.old_snapshot_comparison_used);
        add("Resolved config SHA256", payload.resolved_config_sha256);
        add("Scientific SHA256", payload.scientific_config_sha256);

        holder.appendChild(kv);

        holder.appendChild(el("h3", { text: "Artefacts" }));

        const table = el("table", { class: "compact" }, [
            el("thead", {}, [
                el("tr", {}, [
                    el("th", { text: "Path" }),
                    el("th", { class: "num", text: "Bytes" }),
                    el("th", { text: "SHA256" }),
                ]),
            ]),
        ]);

        const body = el("tbody");

        (release.artifacts || []).forEach(function (artifact) {
            body.appendChild(
                el("tr", {}, [
                    el("td", { class: "mono", text: artifact.path }),
                    el("td", { class: "num", text: num(artifact.bytes) }),
                    el("td", { class: "mono", text: artifact.sha256 }),
                ])
            );
        });

        table.appendChild(body);
        holder.appendChild(table);
    } catch (error) {
        holder.appendChild(
            el("div", { class: "fail-block", text: error.message })
        );
    }
}

/* ----------------------------------------------------------------- runs */

async function loadRuns() {
    const body = $("runs-table").tBodies[0];

    body.innerHTML = "";
    $("run-detail").innerHTML = "";

    try {
        const params = new URLSearchParams();

        if ($("cfg-profile").value) {
            params.set("config", $("cfg-profile").value);
        }

        const payload = await api("/api/runs?" + params.toString());

        if (!payload.runs.length) {
            $("run-detail").appendChild(
                el("p", {
                    class: "note",
                    text: "No runs recorded under " + payload.run_root,
                })
            );
            return;
        }

        payload.runs.forEach(function (run) {
            const row = el("tr", {}, [
                el("td", {}, [
                    el("button", { class: "link", text: run.run_id }),
                ]),
                el("td", { class: "mono", text: run.created_at || "—" }),
                el("td", {}, [stateCell(run.status)]),
                el("td", { class: "mono", text: run.snapshot_id || "—" }),
                el("td", {
                    class: "mono",
                    text: (run.scientific_config_sha256 || "").slice(0, 16),
                }),
                el("td", {
                    class: "mono",
                    text: (run.resolved_config_sha256 || "").slice(0, 16),
                }),
                el("td", { class: "num", text: num(run.stage_count) }),
            ]);

            row.querySelector("button").addEventListener("click", function () {
                loadRunDetail(run.run_id);
            });

            body.appendChild(row);
        });
    } catch (error) {
        $("run-detail").appendChild(
            el("div", { class: "fail-block", text: error.message })
        );
    }
}

async function loadRunDetail(runId) {
    const holder = $("run-detail");

    holder.innerHTML = "";

    try {
        /* Read-only: this endpoint reads the run's own records. It never
         * writes run.json, appends an event, re-resolves a snapshot or
         * launches anything. */
        const payload = await api(
            "/api/runs/" + encodeURIComponent(runId) + "/timeline"
        );

        holder.appendChild(el("h3", { text: "Run " + payload.run_id }));

        holder.appendChild(
            el("p", {
                class: "note",
                text:
                    "Read-only historical record. Stages are shown in "
                    + "canonical scientific order. Select a row to inspect "
                    + "what this run recorded for that stage.",
            })
        );

        /* ---------------------------------------------- run summary */

        const kv = el("dl", { class: "kv" });

        function add(label, value) {
            kv.appendChild(el("dt", { text: label }));
            kv.appendChild(
                el("dd", {
                    text:
                        value === null || value === undefined || value === ""
                            ? "not recorded"
                            : String(value),
                })
            );
        }

        const snapshot = payload.snapshot || {};
        const git = payload.git || {};
        const env = payload.environment || {};

        add("Status", payload.status);
        add("Created", payload.created_at);
        add("Resolved snapshot", snapshot.snapshot_id);
        add("Snapshot selection mode", snapshot.selection_mode);
        add("Scientific SHA256", payload.scientific_config_sha256);
        add("Resolved config SHA256", payload.resolved_config_sha256);
        add("Git branch", git.branch);
        add("Git commit", git.commit);
        add("Working tree dirty", git.working_tree_dirty);
        add("Python", env.python_version);
        add("numpy / scipy / pyarrow",
            [env.numpy_version, env.scipy_version, env.pyarrow_version]
                .filter(Boolean).join("  /  "));
        add("BRI implementation", env.bri_implementation);
        add("Run directory", payload.run_directory);
        add("Configuration file", payload.config_file);
        add("Overrides",
            (payload.config_overrides || []).join(", ") || "none");

        holder.appendChild(kv);

        /* -------------------------------------------- stage timeline */

        holder.appendChild(el("h3", { text: "Pipeline timeline" }));

        const table = el("table", { class: "compact" }, [
            el("thead", {}, [
                el("tr", {}, [
                    el("th", { text: "Canonical stage" }),
                    el("th", { text: "Layer" }),
                    el("th", { text: "Status" }),
                    el("th", { text: "Validation" }),
                    el("th", { class: "num", text: "Input" }),
                    el("th", { class: "num", text: "Output" }),
                    el("th", { text: "Slurm" }),
                    el("th", { text: "Reused" }),
                ]),
            ]),
        ]);

        const body = el("tbody");

        /* payload.timeline is already in canonical order. Never re-sort it. */
        payload.timeline.forEach(function (stage) {
            const nameCell = el("td", {}, [
                el("span", { class: "canonical-label", text: stage.label }),
                document.createTextNode("  "),
                document.createTextNode(stage.title),
            ]);

            if (stage.shared_producer) {
                nameCell.appendChild(el("br"));
                nameCell.appendChild(
                    el("span", {
                        class: "canonical-sub",
                        text: "shares producer " + stage.producer,
                    })
                );
            } else if (stage.parent) {
                nameCell.appendChild(el("br"));
                nameCell.appendChild(
                    el("span", {
                        class: "canonical-sub",
                        text: "part of " + stage.parent,
                    })
                );
            }

            const row = el("tr", {
                class: "stage-row role-" + stage.role,
                "aria-expanded": "false",
                tabindex: "0",
            }, [
                nameCell,
                el("td", {}, [
                    el("span", { class: "layer-tag", text: stage.layer }),
                ]),
                el("td", {}, [stateCell(stage.status)]),
                el("td", {}, [stateCell(stage.validation)]),
                el("td", { class: "num", text: num(stage.input_count) }),
                el("td", { class: "num", text: num(stage.output_count) }),
                el("td", {
                    class: "mono",
                    text: (stage.slurm_job_ids || []).join(", ") || "—",
                }),
                el("td", { text: stage.reused ? "reused" : "—" }),
            ]);

            const detail = el("tr", { class: "stage-detail", hidden: "" }, [
                el("td", { colspan: "8" }, [
                    el("div", { class: "note", text: "loading…" }),
                ]),
            ]);

            function toggle() {
                const open = row.getAttribute("aria-expanded") === "true";

                row.setAttribute("aria-expanded", open ? "false" : "true");
                detail.hidden = open;

                if (!open && !detail.dataset.loaded) {
                    detail.dataset.loaded = "1";
                    loadStageDetail(payload.run_id, stage, detail);
                }
            }

            row.addEventListener("click", toggle);
            row.addEventListener("keydown", function (event) {
                if (event.key === "Enter" || event.key === " ") {
                    event.preventDefault();
                    toggle();
                }
            });

            body.appendChild(row);
            body.appendChild(detail);
        });

        table.appendChild(body);
        holder.appendChild(el("div", { class: "scroll-x" }, [table]));
    } catch (error) {
        holder.appendChild(
            el("div", { class: "fail-block", text: error.message })
        );
    }
}

/* ------------------------------------------ historical stage inspection */

function detailBlock(title, pairs) {
    const dl = el("dl");

    let shown = 0;

    Object.keys(pairs).forEach(function (label) {
        const value = pairs[label];

        if (value === undefined) {
            return;
        }

        shown += 1;

        const missing =
            value === null
            || value === ""
            || value === "not recorded";

        dl.appendChild(el("dt", { text: label }));
        dl.appendChild(
            el("dd", {
                class: missing ? "unrecorded" : null,
                text: missing ? "not recorded" : String(value),
            })
        );
    });

    if (!shown) {
        return null;
    }

    return el("div", { class: "detail-block" }, [
        el("h4", { text: title }),
        dl,
    ]);
}

async function loadStageDetail(runId, stage, container) {
    const cell = container.querySelector("td");

    cell.innerHTML = "";

    try {
        const detail = await api(
            "/api/runs/" + encodeURIComponent(runId)
            + "/stages/" + encodeURIComponent(stage.key)
        );

        const grid = el("div", { class: "detail-grid" });

        const id = detail.identity;
        const st = detail.status;
        const ex = detail.execution;

        [
            detailBlock("Identity", {
                "Canonical stage": id.canonical_label,
                "Canonical name": id.canonical_title,
                "Parent stage": id.parent_stage,
                "Layer": id.layer,
                "Producer": id.producer,
                "Shares producer with":
                    (id.sibling_identities || []).join(", ") || undefined,
                "Implementation": id.implementation,
                "Frozen evidence": id.frozen_output,
            }),
            detailBlock("Status", {
                "Status": st.status,
                "Validation": st.validation,
                "Reused": st.reused ? "yes" : "no",
                "Attempts": st.attempts,
            }),
            detailBlock("Configuration", detail.configuration || {}),
            detailBlock("Inputs", {
                "Input count": (detail.inputs || {}).input_count,
                "Upstream stages":
                    ((detail.inputs || {}).upstream_stages || []).join(", "),
                "Input manifest": (detail.inputs || {}).manifest_path,
            }),
            detailBlock("Outputs", {
                "Output count": (detail.outputs || {}).output_count,
                "Output path": (detail.outputs || {}).output_path,
                "Summary": (detail.outputs || {}).summary_path,
            }),
            detailBlock("Validation", {
                "Verdict": (detail.validation || {}).verdict,
                "Gate": (detail.validation || {}).gate,
            }),
            detailBlock("Execution provenance", {
                "Started": ex.started_at,
                "Finished": ex.finished_at,
                "Runtime (s)": ex.runtime_seconds,
                "Slurm job IDs": (ex.slurm_job_ids || []).join(", "),
                "Entry point": ex.entry_point,
                "Git commit": ex.git_commit,
                "Git branch": ex.git_branch,
                "Working tree dirty": String(ex.working_tree_dirty),
                "Scientific SHA256": ex.scientific_config_sha256,
                "Resolved config SHA256": ex.resolved_config_sha256,
            }),
            detailBlock("Reuse", {
                "Reused": detail.reuse.reused ? "yes" : "no",
                "Why": detail.reuse.explanation,
            }),
        ].forEach(function (block) {
            if (block) {
                grid.appendChild(block);
            }
        });

        cell.appendChild(el("p", { class: "note", text: id.purpose }));

        if (id.note && id.note !== "not recorded") {
            cell.appendChild(el("p", { class: "note", text: id.note }));
        }

        cell.appendChild(grid);

        /* --------------------------------- duplicate / Mol* navigation */

        if (detail.duplicate_navigation) {
            const nav = detail.duplicate_navigation;
            const actions = el("div", { class: "actions" });
            const button = el("button", {
                class: "action",
                text: nav.label,
            });

            button.addEventListener("click", function (event) {
                event.stopPropagation();
                openDuplicateExplorer(nav.filters || {});
            });

            actions.appendChild(button);
            actions.appendChild(
                el("span", { class: "note", text: nav.caveat })
            );
            cell.appendChild(actions);
        }

        /* ------------------------------------------ artefact browser */

        const files = detail.artefacts || [];

        if (files.length) {
            cell.appendChild(
                el("h4", {
                    class: "artefact-list",
                    text: "Output artefacts (" + files.length + ")",
                })
            );

            const list = el("table", { class: "compact" }, [
                el("thead", {}, [
                    el("tr", {}, [
                        el("th", { text: "Artefact" }),
                        el("th", { class: "num", text: "Bytes" }),
                        el("th", { text: "SHA256" }),
                        el("th", { text: "" }),
                    ]),
                ]),
            ]);

            const listBody = el("tbody");

            files.forEach(function (file) {
                const action = el("td");

                if (file.previewable) {
                    const view = el("button", {
                        class: "link",
                        text: "preview",
                    });

                    view.addEventListener("click", function (event) {
                        event.stopPropagation();
                        previewArtefact(file, cell);
                    });

                    action.appendChild(view);
                } else {
                    action.appendChild(
                        el("span", { class: "note", text: "metadata only" })
                    );
                }

                listBody.appendChild(
                    el("tr", {}, [
                        el("td", {
                            class: "mono",
                            text: file.relative_path || file.name,
                        }),
                        el("td", { class: "num", text: num(file.bytes) }),
                        el("td", {
                            class: "mono",
                            text: file.sha256
                                ? file.sha256.slice(0, 16)
                                : "—",
                        }),
                        action,
                    ])
                );
            });

            list.appendChild(listBody);
            cell.appendChild(el("div", { class: "scroll-x" }, [list]));
        }
    } catch (error) {
        cell.appendChild(
            el("div", { class: "fail-block", text: error.message })
        );
    }
}

async function previewArtefact(file, container) {
    const existing = container.querySelector(".artefact-preview");

    if (existing) {
        existing.remove();
    }

    const holder = el("div", { class: "artefact-preview" });

    container.appendChild(holder);
    holder.appendChild(el("div", { class: "note", text: "loading…" }));

    try {
        const payload = await api(
            "/api/artefact?path=" + encodeURIComponent(file.path)
            + "&limit=50"
        );

        holder.innerHTML = "";

        const preview = payload.preview || {};

        if (preview.truncated) {
            holder.appendChild(
                el("div", {
                    class: "preview-note",
                    text:
                        preview.reason
                        || "Preview truncated: showing a bounded sample only.",
                })
            );
        }

        if (preview.kind === "json") {
            holder.appendChild(
                el("pre", {
                    text: JSON.stringify(preview.content, null, 2),
                })
            );
        } else if (preview.kind === "text") {
            holder.appendChild(el("pre", { text: preview.content }));
        } else if (preview.kind === "table" || preview.kind === "parquet") {
            holder.appendChild(renderPreviewTable(preview));
        } else {
            holder.appendChild(
                el("div", {
                    class: "note",
                    text: preview.reason || "No inline preview available.",
                })
            );
        }
    } catch (error) {
        holder.innerHTML = "";
        holder.appendChild(
            el("div", { class: "fail-block", text: error.message })
        );
    }
}

function renderPreviewTable(preview) {
    const columns =
        preview.kind === "parquet"
            ? preview.columns.map(function (c) { return c.name; })
            : preview.columns;

    const head = el("tr", {}, columns.map(function (name, index) {
        const type =
            preview.kind === "parquet"
                ? preview.columns[index].type
                : null;

        return el("th", { title: type, text: name });
    }));

    const rows = (preview.rows || []).map(function (row) {
        const values =
            preview.kind === "parquet"
                ? columns.map(function (name) { return row[name]; })
                : row;

        return el("tr", {}, values.map(function (value) {
            return el("td", {
                class: "mono",
                text:
                    value === null || value === undefined
                        ? "—"
                        : String(value),
            });
        }));
    });

    const wrapper = el("div", { class: "scroll-x" }, [
        el("table", { class: "compact" }, [
            el("thead", {}, [head]),
            el("tbody", {}, rows),
        ]),
    ]);

    if (preview.kind === "parquet") {
        const meta = el("div", { class: "note" });

        meta.textContent =
            num(preview.row_count) + " rows, "
            + preview.columns.length + " columns; showing "
            + preview.row_preview_count + ".";

        wrapper.insertBefore(meta, wrapper.firstChild);
    }

    return wrapper;
}

/* Jump to the Duplicate Explorer with the stage's filters applied. Mol* is
 * reached from there; neither ever reclassifies a pair. */
function openDuplicateExplorer(filters) {
    if (filters.relationship) {
        $("dup-rel").value = filters.relationship;
    }

    if (filters.exact_only) {
        $("dup-class").value = "exact";
    }

    showView("duplicates");
    searchDuplicates(0);
}

/* ------------------------------------------------------------ start run */

async function startRun() {
    const holder = $("run-result");

    holder.innerHTML = "";

    try {
        const payload = await api("/api/run", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify(requestBody()),
        });

        holder.appendChild(
            el("div", { class: "warn-block" }, [
                el("div", { text: "Run created: " + payload.run_id }),
                el("div", {
                    class: "note",
                    text:
                        "Provenance was written before any work. Long stages "
                        + "are submitted from the CLI on a login node, not "
                        + "from the browser.",
                }),
            ])
        );

        const kv = el("dl", { class: "kv" }, [
            el("dt", { text: "Run directory" }),
            el("dd", { text: payload.run_directory }),
            el("dt", { text: "Resolved config SHA256" }),
            el("dd", { text: payload.resolved_config_sha256 }),
            el("dt", { text: "CLI equivalent" }),
            el("dd", { text: payload.cli_equivalent }),
        ]);

        holder.appendChild(kv);

        if (payload.commands.length) {
            holder.appendChild(el("h3", { text: "Outstanding stage commands" }));

            payload.commands.forEach(function (entry) {
                holder.appendChild(
                    el("pre", {
                        text:
                            "# " + entry.stage_id + "\n"
                            + (entry.error
                                ? "ERROR: " + entry.error
                                : entry.argv.join(" \\\n    ")),
                    })
                );
            });
        } else {
            holder.appendChild(
                el("p", {
                    class: "note",
                    text:
                        "Every stage already has validated output for this "
                        + "configuration; nothing needs to run.",
                })
            );
        }

        toast("Run " + payload.run_id + " created.");
    } catch (error) {
        holder.appendChild(
            el("div", { class: "fail-block", text: error.message })
        );
        toast(error.message, true);
    }
}

/* ---------------------------------------------------------------- about */

function renderAbout() {
    const body = $("about-stages").tBodies[0];

    body.innerHTML = "";

    /* The canonical timeline, in canonical order. Never re-sorted. */
    state.bootstrap.canonical_timeline.forEach(function (stage) {
        const name = el("td", { class: "canonical" }, [
            el("span", { class: "canonical-label", text: stage.label }),
            document.createTextNode("  " + stage.title),
        ]);

        const annotations = [];

        if (stage.shared_producer) {
            annotations.push("shares producer with a sibling stage");
        }

        if (stage.parent) {
            annotations.push("part of " + stage.parent);
        }

        if (!stage.producer) {
            annotations.push("not executed by the orchestrated pipeline");
        }

        if (annotations.length) {
            name.appendChild(el("br"));
            name.appendChild(
                el("span", {
                    class: "canonical-sub",
                    text: annotations.join(" · "),
                })
            );
        }

        const purpose = el("td", { text: stage.purpose });

        if (stage.note) {
            purpose.appendChild(el("br"));
            purpose.appendChild(
                el("span", { class: "canonical-sub", text: stage.note })
            );
        }

        body.appendChild(
            el("tr", { class: "role-" + stage.role }, [
                name,
                el("td", {}, [
                    el("span", { class: "layer-tag", text: stage.layer }),
                ]),
                el("td", {
                    class: "mono",
                    text: stage.producer || "—",
                }),
                purpose,
            ])
        );
    });
}

/* ------------------------------------------------------------------ init */

async function init() {
    state.bootstrap = await api("/api/bootstrap");

    const profile = $("cfg-profile");

    profile.appendChild(
        el("option", {
            value: "",
            text: "built-in validated defaults (no file)",
        })
    );

    state.bootstrap.profiles.forEach(function (entry) {
        profile.appendChild(
            el("option", { value: entry.path, text: entry.path })
        );
    });

    if (state.bootstrap.default_config_path) {
        profile.value = state.bootstrap.default_config_path;
    }

    fillFormFromDefaults(state.bootstrap.defaults);
    renderAbout();

    document.querySelectorAll("nav button").forEach(function (button) {
        button.addEventListener("click", function () {
            showView(button.dataset.view);
        });
    });

    $("btn-resolve").addEventListener("click", function () {
        resolveConfiguration(false);
    });

    $("btn-reset").addEventListener("click", function () {
        fillFormFromDefaults(state.bootstrap.defaults);
        $("cfg-snapshot").value = "";
        $("cfg-snapshot-manual").value = "";
        resolveConfiguration(false);
    });

    $("btn-plan-from-config").addEventListener("click", async function () {
        await resolveConfiguration(true);
        showView("pipeline");
        loadPlan();
    });

    $("btn-list-snapshots").addEventListener("click", listSnapshots);
    $("btn-start-run").addEventListener("click", startRun);

    $("btn-dup-search").addEventListener("click", function () {
        searchDuplicates(0);
    });

    $("btn-dup-clear").addEventListener("click", function () {
        ["dup-pdb", "dup-chain", "dup-minlen", "dup-maxlen", "dup-mind",
            "dup-maxd"].forEach(function (id) {
            $(id).value = "";
        });

        $("dup-class").value = "";
        $("dup-rel").value = "";

        searchDuplicates(0);
    });

    $("dup-prev").addEventListener("click", function () {
        const limit = Number($("dup-limit").value);

        searchDuplicates(Math.max(0, state.dupOffset - limit));
    });

    $("dup-next").addEventListener("click", function () {
        const limit = Number($("dup-limit").value);

        searchDuplicates(state.dupOffset + limit);
    });

    initTheme();

    await resolveConfiguration(true);
}

init().catch(function (error) {
    document.body.innerHTML =
        "<pre style=\"margin:16px\">UI failed to start: "
        + String(error.message)
        + "</pre>";
});
