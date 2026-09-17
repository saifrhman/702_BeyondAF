# Sequence-redundancy resolution (Stage 15)

What the stage does, why it exists, why its parameters have the values they do,
and exactly how many chains came out.

Every figure below is traceable to a file in this repository or to a recorded
stage run; the sources are named at the end of each section.

---

## 1. Why the stage exists

The PDBClean pipeline removes chains that are redundant **by shape**. It says
nothing about whether two chains are redundant **by sequence**, and the curated
population was never reduced on that axis: MMseqs2 appears elsewhere in this
project only as an alignment *search* tool, generating one MSA per retained
chain against UniRef30. It was never used to cluster the training set against
itself.

So the geometric release

    PDBClean-20260101-protocol3.2-comp702-v1-dedup-v1      499,770 chains

contains substantial sequence redundancy by design. Those 499,770 chains carry
only **142,056 distinct sequences**: on average each sequence appears about 3.5
times. Training on that population shows the model the same sequence many times
over, weighted by how often crystallographers happened to deposit it.

Stage 15 asks the second question, after geometry has answered the first.

Sources: `outputs/releases/PDBClean-20260101-protocol3.2-comp702-v1-dedup-v1/data/retained_chains.parquet`;
distinct-sequence count from the Stage-15 run summary.

---

## 2. What the stage actually does

1. Takes the retained chains of a completed Gold release, with their
   **retained (post-trimming) sequences** — the same sequences the alignment
   stage used, not the raw deposited ones.
2. Clusters them with MMseqs2.
3. Keeps **one chain per cluster**, chosen by this project's own ranking.
4. Publishes the survivors as a new, separately identified release. The
   geometry-only release it consumes is read-only input and is never modified.

### Representative selection is not delegated to MMseqs2

This is the part that matters most, and it is deliberate.

MMseqs2's `--cluster-mode 0` is a greedy set cover. It is order-dependent, and
on this very population it was observed placing **byte-identical sequences under
different representatives** — `2n0k_A` and `2n0k_B`, the same 89-residue
sequence, landed under `2wj7_E` and `4jut_C` respectively.

MMseqs2 therefore decides cluster **membership only**. The survivor within each
cluster is chosen by the deterministic ranking Stage 14b already applies to
geometric components:

1. untrimmed preferred over terminal-trimmed
2. fewer defective (dirty) residues preferred
3. better nominal resolution preferred — but only where the cluster's
   experimental methods are comparable
4. canonical chain key, as a total deterministic tie-break

The ranking helpers are **imported from the Stage 14b entry point**
(`scripts/select_stage14_representatives.py`) rather than copied, so the two
stages cannot drift apart in what they consider a good representative.

**Caveat on rule 3.** The resolution term only fires when every member of a
cluster shares one experimental-method tuple, and deposition metadata exists for
only **6,357 of 118,197 entries (5.4%)** — the `downstream_metadata` stage
fetches only the depositions that participate in *geometric* near-duplicate
pairs. In practice rule 3 applied to **4,717 of 84,156 clusters (5.6%)**. For the
remaining 94% selection falls through to rules 1, 2 and 4, which are defined for
every chain and are themselves deterministic. Selection is therefore
reproducible, but it is not resolution-aware at population scale, and the
dissertation should not claim that it is.

Sources: `src/pdbclean/sequence_clustering.py`;
`outputs/pdbclean/20260101/protocol3.2-comp702-v1/sequence_clustering/global_summary.json`.

---

## 3. The parameters, and why

    tool          MMseqs2 18.8cc5c
    subcommand    easy-cluster
    min_seq_id    1.0
    coverage      0.8
    cov_mode      0      (bidirectional)
    cluster_mode  0      (set cover)

### `min_seq_id 1.0`

At 1.0 the stage removes only chains whose sequence is indistinguishable from
their representative over the covered region. This is the redundancy nobody
disputes. Every looser threshold starts discarding genuine homologues that a
structure-prediction model has reason to see — and a homologue is training
signal, not noise.

It is the conservative default for a training population. The value is
configuration, not a constant, so a threshold study is a one-key change.

### `coverage 0.8`, `cov_mode 0`

`cov_mode 0` requires the alignment to cover at least 80% of **both** sequences.
This is the parameter that stops a short chain being absorbed into a long one on
the strength of a shared domain: under `cov_mode 0` a 90-residue chain cannot
cluster with a 400-residue chain regardless of identity. One-sided coverage
modes (1, 2) would inflate the removals with partial matches that are not the
same molecule.

### `easy-cluster`, not `linclust`

`linclust` is linear-time and much faster, but finds far fewer remote
homologies — exactly in the 30–50% identity range a threshold sweep turns on.
Using it would understate the "same cluster" outcome and so *overstate* how much
geometric de-duplication is orthogonal to sequence clustering. That is the
direction of error to avoid here, so the slower, more sensitive mode is used.

Sources: `src/pdbclean/defaults.py`,
`config/pdbclean/profiles/comp702_seqclust_20260101.yaml`.

---

## 4. What the threshold sweep showed

Eight settings, each varying exactly one parameter on the same 499,770-chain
population, same clusterer, same ranking.

    min_seq_id (coverage 0.8)      retained    removed   singletons   largest
      1.00                           84,156    415,614      35,783     2,537
      0.90                           62,792    436,978      25,618     3,119
      0.70                           52,610    447,160      20,993     5,055
      0.50                           44,799    454,971      17,917     5,741
      0.30                           36,322    463,448      14,896     7,097

    coverage (min_seq_id 1.0)      retained    removed   singletons   largest
      0.80                           84,156    415,614      35,783     2,537
      0.90                           87,255    412,515      38,029     2,537
      0.95                           92,319    407,451      41,866     2,537
      1.00                          143,226    356,544      84,579     2,537

### The finding that matters

The population contains **142,056 distinct sequences**, yet clustering at
`min_seq_id 1.0, -c 0.8` returns only **84,156** clusters. Roughly 58,000
sequences were merged with something they are not byte-identical to.

The coverage sweep explains it. At `-c 1.0` the retained count is **143,226** —
essentially recovering exact-sequence de-duplication. The gap at `-c 0.8` is
**length-variant collapse**: two chains with the same core but different
observed termini satisfy 80% bidirectional coverage and merge.

The identity distributions corroborate this. At `cov 1.0` every one of the
356,544 removals is 100% identical to its representative. At `cov 0.8` there are
**513 removals below 30% identity** — pairs that are not the same sequence at
all, merged on an 80%-covered alignment.

### Consequence, stated plainly

**The dataset now training uses `coverage 0.8`, and the sweep evidence argues
`coverage 1.0` is the more defensible choice.** At `min_seq_id 1.0, -c 1.0` the
operation is exactly "remove chains whose sequence is byte-identical to their
representative over the full length" — no length variants collapsed, no
sub-30%-identity removals to explain in a viva.

The training population described in §6 is the `-c 0.8` one. If the argument in
the dissertation is "we removed exact sequence duplicates", `-c 1.0` matches
that claim and `-c 0.8` does not. This is a live decision, not a settled one.

Source: `COMP702_followup/sweep/*.json` (eight recorded runs).

---

## 5. Validation

Gates enforced by the stage, all passing:

    chain accounting reconciles exactly           84,156 + 415,614 = 499,770
    every removed chain resolves to exactly one representative
    that representative is in the retained set
    no representative is itself removed
    every retained chain is its own representative

**Determinism was measured, not asserted.** The stage was run twice as
independent Slurm jobs (`10561416`, `10561769`) on the same input:

    finalized/retained_chains.parquet       e0a4b0fea7a64ced0c7be0d5   IDENTICAL
    finalized/removed_chain_audit.parquet   07944c74cb4a30028f11b527   IDENTICAL
    global_summary.json                                                identical

This is not luck: cluster ids and their members are sorted before selection, so
the result cannot depend on the order MMseqs2 emitted them in.

### Identity of each removal

Every removed chain records its pairwise identity to its representative, and the
method used to obtain it:

    mmseqs_alignment       404,415     alignment identity from MMseqs2
    exact_match              8,460     byte-identical, 1.0 by inspection
    ungapped_best_offset     2,528     ungapped scan, a LOWER BOUND
    ungapped_positional        211     equal length, position by position

    identity  100%   414,305  (99.685%)
              >=90%      771
              >=70%       20
              <50%       518

MMseqs2's k-mer prefilter cannot seed very short sequences, so 11,199 pairs came
back with no alignment (median length 10 residues). Rather than leave those
null, they are resolved exactly where possible. **The 518 pairs below 50%
identity all come from the two fallback methods, never from an MMseqs2
alignment** — an ungapped comparison underestimates whenever the true alignment
has gaps, so treat those as a floor, not a measurement.

Sources: `sequence_clustering/global_summary.json`;
`sequence_clustering/finalized/removed_chain_audit.parquet`
(column `seqdedup_identity_method`).

---

## 6. The dataset being trained

    input   (geometric release)                    499,770 chains
    removed (sequence-redundant)                   415,614 chains
    RETAINED — the training population               84,156 chains

    distinct sequences among them                    83,140
    distinct PDB entries                             (subset of 118,197)
    clusters formed                                  84,156
    largest cluster                                   2,537 chains
    cluster sizes   1: 35,783 | 2: 17,328 | 3-10: 24,113
                    11-100: 6,389 | 101-1000: 537 | >1000: 6

Release identifier:

    PDBClean-20260101-protocol3.2-comp702-v1-dedup-v1-seqid100-v1

That is an **83.2% reduction** from the geometric population. It is a large
intervention, and the comparison against the geometry-only run will be
confounded by dataset size as much as by redundancy — worth stating in any
comparison of the two models.

### Training-ready inputs

MSAs are content-addressed by sequence SHA256, and the survivors are a subset of
a population whose MSAs already exist, so alignment generation was a
re-projection rather than a regeneration. This was verified, not assumed:

    distinct sequences needing an MSA                83,140
    MSAs already present                             83,140   (100.000%)
    MSAs generated                                        0

Coverage gate (`verify_alignment_index.py`, 2,000 chains sampled):

    chains in index                                  84,156
    chains without an index entry                         0
    entries with no chain                                 0
    unreadable entries                                    0
    wrong sequence delivered                              0
    STATUS: PASS

Sources:
`outputs/openfold_training/20260101/pdbclean-dedup-v1-seqid100-v1/training_manifest.json`,
`alignment_index_verification.json`.

---

## 7. Known limitations

**The result is not exactly sequence-deduplicated.** The retained population
still contains **1,016 chains whose sequence duplicates another retained chain**,
across 105 sequences (median length 7 residues, longest 130 — for example
`4v9h_EB`, `4v9i_IA` and `4v9i_OC` all carry the same 130-residue sequence).
This is the MMseqs2 non-transitivity described in §2 surviving into the output.
An exact post-pass over byte-identical sequences would close it; that would
change the population, so it was not applied without asking.

**The resolution tie-break is largely inert** — see §2. Selection is
deterministic but not resolution-aware for 94% of clusters.

**518 identity values are lower bounds, not alignments** — see §5.

**Coverage 0.8 collapses length variants** — see §4. This is the most
consequential open question about the dataset.

---

## 8. Reproducing it

```bash
pdbclean run --config config/pdbclean/profiles/comp702_seqclust_20260101.yaml
```

or the stage alone:

```bash
python -m pdbclean.sequence_clustering_production \
    --config <frozen stage config> --pipeline-git-commit <sha>
```

The stage is registered as **Stage 15**, `depends_on=("gold_release",)`, and
appears in `pdbclean stages`, `pdbclean plan` and `pdbclean run` like every
other stage. It is **disabled by default**: a run that does not ask for it
publishes the geometry-only population, and the planner reports the stage as
*not applicable* rather than as work forever outstanding. That switch,
`sequence_clustering.enabled`, is the only difference between producing the two
datasets.

To project the result into OpenFold inputs:

```bash
python scripts/openfold_training/build_seqclust_training_view.py \
    --release  outputs/releases/PDBClean-...-dedup-v1-seqid100-v1 \
    --msa-store <msa store> \
    --output-root outputs/openfold_training/20260101/pdbclean-dedup-v1-seqid100-v1
```
