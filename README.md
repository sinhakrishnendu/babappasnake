# BABAPPASnake
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19048331.svg)](https://doi.org/10.5281/zenodo.19048331)

`babappasnake` is a command-line workflow for episodic positive selection analysis on one orthogroup at a time. It is built for practical comparative genomics: reproducible run directories, guided stepwise execution, resumable restarts after interruption, and robustness summaries across multiple alignment and trimming pathways.

This README is the user manual for the current CLI and workflow behavior.

## Contents

1. [What The Pipeline Does](#what-the-pipeline-does)
2. [Installation And Environment](#installation-and-environment)
3. [Input Requirements](#input-requirements)
4. [Orthogroup Discovery Modes](#orthogroup-discovery-modes)
5. [Running The Pipeline](#running-the-pipeline)
6. [Resume And Recovery](#resume-and-recovery)
7. [Workflow Stages](#workflow-stages)
8. [Alignment And Robustness Design](#alignment-and-robustness-design)
9. [Recombination, Selection, And ASR](#recombination-selection-and-asr)
10. [Outputs And Run Directory Layout](#outputs-and-run-directory-layout)
11. [CLI Reference](#cli-reference)
12. [Troubleshooting](#troubleshooting)
13. [Developer Notes](#developer-notes)
14. [License](#license)

## What The Pipeline Does

A typical full run performs these stages:

1. Define a one-to-one orthogroup from a protein query and a proteome panel.
2. Map user-provided CDS to the selected orthogroup proteins with coding-quality filters.
3. Build protein and codon alignments across one or more alignment methods.
4. Run both `raw` and `clipkit` pathway variants for robustness.
5. Infer pathway-specific trees with IQ-TREE.
6. Optionally root trees with an outgroup text query.
7. Optionally screen for recombination with HyPhy GARD.
8. Run HyPhy aBSREL and MEME.
9. Select branch-site foregrounds dynamically from aBSREL output.
10. Run branch-site codeml.
11. Optionally run pathway-level codeml ASR.
12. Optionally extract ancestor and descendant sequences plus branch substitutions.
13. Write pathway summaries, robustness reports, and run provenance.

Design assumptions:

- The workflow operates on one orthogroup at a time.
- CDS can be unavailable at the start. The pipeline supports a protein-first checkpointed workflow.
- Outgroup is optional. If it is absent or unusable, downstream analysis continues with the unrooted tree.
- Robustness mode is always enforced internally as `raw + clipkit`, even if you request a single trimming mode.

## Installation And Environment

### Recommended one-command environment

```bash
conda create -n babappasnake -c conda-forge -c bioconda \
  python=3.11 blast orthofinder iqtree hyphy paml clipkit mafft prank pip
conda activate babappasnake
pip install babappasnake
```

### External tools by feature

Always needed for a complete end-to-end run:

- `blastp`
- `makeblastdb`
- `iqtree` or `iqtree2` or `iqtree3`
- `hyphy`
- `codeml` from `paml`
- `clipkit`

Needed only for selected orthogroup modes:

- `orthofinder` for `--orthogroup-method orthofinder`
- `orthofinder` for `--orthogroup-method rbh_fallback`

Needed only for selected alignment methods:

- `mafft` when `--alignment-methods` is `2` or `4`
- `prank` when `--alignment-methods` is `3` or `4`

Python-side note:

- `babappalign` is installed automatically as a dependency of `babappasnake`.

### Verify installation

```bash
babappasnake --help
which blastp makeblastdb iqtree iqtree2 iqtree3 hyphy codeml clipkit
which orthofinder mafft prank
which babappasnake
```

## Input Requirements

### `--prot`

Directory of protein FASTA files, one file per species.

Requirements and behavior:

- Supported extensions are `.fa`, `.faa`, and `.fasta`, case-insensitive.
- Hidden files and macOS metadata sidecars such as `._*` and `.DS_Store` are ignored.
- Empty or malformed FASTA files are skipped with warnings.

### `--query`

Protein FASTA containing exactly one query sequence.

### `--cds`

Optional on the first run.

Behavior:

- If CDS is supplied at the start, the workflow can proceed through codon, tree, HyPhy, codeml, and summary stages in one run.
- If CDS is not supplied, the workflow stops after orthogroup definition and writes `orthogroup/WAITING_FOR_CDS.txt`.
- You can then place CDS at `OUTDIR/user_supplied/orthogroup_cds.fasta` and resume.

### `--outgroup`

Optional text query used to root pathway trees by case-insensitive substring matching against tip labels.

Behavior:

- If omitted, rooting is skipped and the unrooted tree is propagated downstream.
- If provided but unmatched or too broad, rooting falls back safely to the unrooted tree instead of crashing the run.

## Orthogroup Discovery Modes

### Default: `--orthogroup-method rbh`

RBH-only mode.

Behavior:

1. Run reciprocal best-hit ortholog recovery.
2. Retain strict one-to-one orthologs only.
3. Stop with an explicit error if no usable strict one-to-one orthogroup is recovered.

This is the default and recommended starting mode.

### Direct OrthoFinder: `--orthogroup-method orthofinder`

OrthoFinder-only mode.

Behavior:

1. Run OrthoFinder.
2. BLAST the query against OrthoFinder orthogroup-member proteins.
3. Rank orthogroups by query support.
4. Retain strict one-to-one orthologs only.

### Comparison fallback: `--orthogroup-method rbh_fallback`

This preserves the older compare-and-choose behavior.

Behavior:

1. Run RBH.
2. Run OrthoFinder.
3. Compute strict one-to-one support for both results.
4. Keep the result with the stronger strict one-to-one count.
5. If tied, keep RBH deterministically.

### OrthoFinder query mapping details

The workflow does not assume that the original query ID is directly present inside an OrthoFinder orthogroup. Instead it does a query-to-members mapping step:

1. Parse `Orthogroups.tsv`.
2. Build a combined FASTA of orthogroup-member proteins.
3. BLAST the query against that combined FASTA.
4. Filter by query coverage threshold.
5. Rank orthogroups by cross-species support and bitscore.
6. Extract the top supported orthogroup and apply strict one-to-one filtering.

### Strict one-to-one rule

A species contributes to the final orthogroup only when exactly one ortholog is retained for that species. This prevents duplicate same-species entries in downstream alignments and trees.

## Running The Pipeline

### Guided interactive mode

```bash
babappasnake
```

Behavior:

- Prompts for run settings.
- Executes the workflow step by step.
- Shows step descriptions, previews, and output summaries.
- Lets you choose `run`, `skip`, or `stop` at each step where appropriate.

Important defaults in guided mode:

- Orthogroup method defaults to `rbh`.
- Alignment methods default to `4` which means all three aligners.
- Trimming is forced internally to robustness mode: `raw + clipkit`.
- Recombination defaults to `none`.
- ASR defaults to `yes`.
- CDS is requested only after orthogroup definition.
- Outgroup is requested after CDS handling and remains optional.

### Non-interactive batch mode

Use this for scripted or HPC-style runs:

```bash
babappasnake \
  --prot /path/to/proteomes \
  --query /path/to/query.fasta \
  --cds /path/to/orthogroup_cds.fasta \
  --outdir run01 \
  --threads 12 \
  --interactive no \
  --guided no
```

### Example: RBH-only run with ASR disabled

```bash
babappasnake \
  --prot /path/to/proteomes \
  --query /path/to/query.fasta \
  --cds /path/to/orthogroup_cds.fasta \
  --orthogroup-method rbh \
  --run-asr no \
  --outdir run_rbh_no_asr \
  --threads 12 \
  --interactive no \
  --guided no
```

### Example: OrthoFinder-only run

```bash
babappasnake \
  --prot /path/to/proteomes \
  --query /path/to/query.fasta \
  --cds /path/to/orthogroup_cds.fasta \
  --orthogroup-method orthofinder \
  --outdir run_orthofinder \
  --threads 12 \
  --interactive no \
  --guided no
```

### Example: enable GARD

```bash
babappasnake \
  --prot /path/to/proteomes \
  --query /path/to/query.fasta \
  --cds /path/to/orthogroup_cds.fasta \
  --recombination gard \
  --gard-mode Faster \
  --gard-rate-classes 3 \
  --outdir run_gard \
  --threads 12 \
  --interactive no \
  --guided no
```

### Two-stage run when CDS is not available initially

Stage 1:

```bash
babappasnake \
  --prot /path/to/proteomes \
  --query /path/to/query.fasta \
  --outdir run01 \
  --interactive no \
  --guided no
```

After orthogroup definition, provide the CDS file and continue:

```bash
babappasnake --resume --outdir run01 --cds /path/to/orthogroup_cds.fasta
```

## Resume And Recovery

Use resume whenever a run stops unexpectedly or intentionally:

```bash
babappasnake --resume --outdir run01
```

Typical reasons to resume:

- workflow failure
- external tool failure
- terminal closure
- power interruption
- machine restart
- manual stop after the orthogroup checkpoint

What `--resume` does:

- reloads `OUTDIR/config.yaml`
- reuses saved analysis settings
- clears stale Snakemake lock state for that run directory
- continues incomplete work instead of restarting from scratch

Guided resume behavior:

- saved and detected completed steps are skipped automatically
- the session restarts near the interruption point

Non-guided resume behavior:

- Snakemake reruns incomplete work only

Files used for resume:

- `OUTDIR/config.yaml`
- `OUTDIR/.babappasnake/resume_state.json`

Allowed resume-time overrides:

- `--cds PATH`
- `--outgroup TEXT`
- `--threads INT`
- `--guided {yes,no}`
- `--snake-args "..." `

Analysis settings that define workflow structure are intentionally not changeable on resume. For example, do not expect `--resume` to accept a new orthogroup method, new alignment mode, new recombination mode, or new ASR setting for an existing run directory.

## Workflow Stages

1. `rbh_orthogroup`
   Build the selected orthogroup by the configured orthogroup backend.
2. `map_cds`
   Map CDS to selected proteins and filter low-quality CDS.
3. `align_proteins_all_methods`
   Build protein alignments.
4. `align_cds_all_methods`
   Build codon alignments.
5. `prepare_branch_inputs_all_pathways`
   Expand each method into `raw` and `clipkit` analysis pathways.
6. `gard_all_pathways`
   Optional recombination screening.
7. `iqtree_ml_all_pathways`
   Infer pathway trees.
8. `root_iqtree_outgroup_all_pathways`
   Root trees when possible, otherwise propagate unrooted trees.
9. `hyphy_exploratory_all_pathways`
   Run aBSREL and MEME.
10. `parse_foregrounds_all_pathways`
    Select significant foreground branches from aBSREL.
11. `prepare_foreground_trees_all_pathways`
    Build branch-labeled trees for branch-site codeml.
12. `branchsite_batch_all_pathways`
    Run branch-site codeml and BH-correct results.
13. `codeml_asr_all_pathways`
    Optional pathway-level ASR.
14. `extract_selected_branch_ancestors`
    Optional branch-level ancestral and descendant sequence extraction.
15. `final_summary_all_pathways`
    Write pathway-specific episodic selection summaries.
16. `robustness_reports`
    Write cross-pathway robustness outputs.
17. `write_run_provenance`
    Write machine-readable provenance for the run.

## Alignment And Robustness Design

### Alignment methods

`--alignment-methods` choices:

- `1`: `babappalign`
- `2`: `mafft`
- `3`: `prank`
- `4`: `babappalign`, `mafft`, and `prank`

### Trimming behavior

The workflow accepts `--trim-strategy`, but runtime robustness mode is always enforced internally as:

- `raw`
- `clipkit`

This means each selected alignment method expands into two pathways.

Example for `--alignment-methods 4`:

- `babappalign_raw`
- `babappalign_clipkit`
- `mafft_raw`
- `mafft_clipkit`
- `prank_raw`
- `prank_clipkit`

This is deliberate. Robustness outputs assume both untrimmed and trimmed pathway variants are available.

## Recombination, Selection, And ASR

### Tree inference and rooting

- IQ-TREE runs once per active pathway.
- Rooting is optional.
- If outgroup matching fails, the workflow continues with the unrooted tree.

### Optional recombination screening with GARD

Controlled by `--recombination {none,gard,auto}`:

- `none`: do not run GARD
- `gard`: run HyPhy GARD
- `auto`: currently an alias of `gard`

Additional controls:

- `--gard-mode {Normal,Faster}`
- `--gard-rate-classes INT`

GARD outputs:

- `recombination/<method>/<trim_state>/gard/gard.json`
- `recombination/<method>/<trim_state>/gard/gard_summary.json`
- `recombination/<method>/<trim_state>/gard/gard.stdout.txt`
- `recombination/<method>/<trim_state>/gard/gard.stderr.txt`

Current interpretation:

- GARD is used as a screening and reporting layer.
- Downstream branch-site codeml is still run on the full-length pathway alignment unless fragment-aware routing is added in a future release.

### HyPhy and branch-site foreground selection

- aBSREL and MEME run per pathway.
- Foregrounds are selected dynamically from aBSREL.
- Default aBSREL dynamic threshold settings:
  - start `0.05`
  - increment `0.01`
  - cap `0.2`
- Selected foregrounds are passed to branch-site codeml.
- Branch-site results are BH-corrected.

### Optional ASR block

When `--run-asr yes`, the workflow:

1. runs codeml ASR per pathway
2. maps selected branches to parent and child nodes
3. extracts ancestor and descendant CDS and amino-acid sequences
4. computes branch substitutions
5. annotates overlaps with MEME and BEB where available

When `--run-asr no`:

- the ASR extraction block is skipped entirely
- HyPhy, branch-site, summaries, robustness reports, and run provenance still complete
- per-pathway `asr_done.json` files are still written with `status: skipped`
- the robustness matrix reports ASR as not completed rather than falsely completed

## Outputs And Run Directory Layout

All outputs live inside `--outdir`.

### Early-stage outputs

- `inputs/query.fasta`
- `inputs/proteomes/`
- `orthogroup/orthogroup_proteins.fasta`
- `orthogroup/orthogroup_headers.txt`
- `orthogroup/rbh_summary.tsv`
- `orthogroup/WAITING_FOR_CDS.txt`
- `user_supplied/orthogroup_cds.fasta`
- `mapped_cds/cds_protein_mapping.tsv`

### Alignment and tree outputs

- `alignments/<method>/orthogroup_proteins.protein.aln.fasta`
- `alignments/<method>/mapped_orthogroup_cds.codon.aln.fasta`
- `alignments/<method>/<trim_state>/orthogroup_proteins.analysis.fasta`
- `alignments/<method>/<trim_state>/mapped_orthogroup_cds.analysis.fasta`
- `tree/<method>/<trim_state>/orthogroup.treefile`
- `tree/<method>/<trim_state>/orthogroup.rooted.treefile`

### Selection outputs

- `hyphy/<method>/<trim_state>/absrel.json`
- `hyphy/<method>/<trim_state>/meme.json`
- `hyphy/<method>/<trim_state>/hyphy_done.json`
- `hyphy/<method>/<trim_state>/significant_foregrounds.tsv`
- `hyphy/<method>/<trim_state>/foreground_threshold.json`
- `branchsite/<method>/<trim_state>/foreground_trees.tsv`
- `branchsite/<method>/<trim_state>/branchsite_results.tsv`

### ASR outputs

Always present per pathway:

- `asr/<method>/<trim_state>/asr_done.json`

Produced only when `--run-asr yes`:

- `asr/<method>/<trim_state>/mlc_asr.txt`
- `asr/<method>/<trim_state>/rst`
- `asr/branch_to_nodes.tsv`
- `asr/ancestor_sequences_cds.fasta`
- `asr/ancestor_sequences_aa.fasta`
- `asr/descendant_sequences_cds.fasta`
- `asr/descendant_sequences_aa.fasta`
- `asr/branch_substitutions.tsv`
- `asr/selected_branch_asr_summary.tsv`
- `asr/asr_extraction_provenance.json`
- `asr/asr_done.json`

### Summary and provenance outputs

- `summary/<method>/<trim_state>/episodic_selection_summary.txt`
- `summary/robustness_matrix.tsv`
- `summary/robustness_consensus.tsv`
- `summary/robustness_narrative.txt`
- `summary/comparative_reproducibility_summary.txt`
- `summary/robustness_publication_table.tex`
- `summary/run_provenance.json`
- `summary/episodic_selection_summary.txt`

### Internal resume state

- `.babappasnake/resume_state.json`
- `.snakemake/`

## CLI Reference

Basic form for a fresh run:

```text
babappasnake --prot PROTEOMES_DIR --query QUERY_FASTA [options]
```

Basic form for a resumed run:

```text
babappasnake --resume --outdir RUN_DIR [resume overrides]
```

### Core options

- `--prot PATH`
- `--query PATH`
- `--cds PATH`
- `--outdir PATH`
- `--interactive {yes,no}`
- `--guided {yes,no}`
- `--resume`
- `--snake-args "..."`

### Orthogroup and alignment options

- `--orthogroup-method {rbh,orthofinder,rbh_fallback}`
- `--coverage FLOAT`
- `--alignment-methods {1,2,3,4}`
- `--trim-strategy {raw,clipkit,both}`
- `--use-clipkit {yes,no}`
- `--clipkit-mode-protein TEXT`
- `--clipkit-mode-codon TEXT`

Practical note:

- `--trim-strategy` and `--use-clipkit` are accepted, but runtime robustness mode is still forced to `both`.

### Tree and selection options

- `--outgroup TEXT`
- `--threads INT`
- `--iqtree-bootstrap INT`
- `--iqtree-bnni {yes,no}`
- `--iqtree-model TEXT`
- `--recombination {none,gard,auto}`
- `--gard-mode {Normal,Faster}`
- `--gard-rate-classes INT`
- `--absrel-branches TEXT`
- `--meme-branches TEXT`
- `--codeml-codonfreq INT`
- `--run-asr {yes,no}`
- `--absrel-p FLOAT`
- `--absrel-dynamic-start FLOAT`
- `--absrel-dynamic-step FLOAT`
- `--absrel-dynamic-max FLOAT`
- `--meme-p FLOAT`

## Troubleshooting

### Run stops at `WAITING_FOR_CDS.txt`

This is expected when CDS was not supplied yet.

Fix:

1. provide `OUTDIR/user_supplied/orthogroup_cds.fasta`
2. run `babappasnake --resume --outdir OUTDIR`

### Outgroup was not available at the start

This is supported.

Behavior:

- the run does not need to crash
- downstream continues with the unrooted tree when no usable outgroup is present

### Interrupted run after crash, reboot, or power loss

Use:

```bash
babappasnake --resume --outdir OUTDIR
```

The workflow will reload saved settings, clear stale lock state, and continue incomplete work.

### OrthoFinder or RBH finds no usable orthogroup

The pipeline stops explicitly when strict one-to-one recovery fails.

Check:

- query quality
- proteome quality
- taxon sampling
- whether the query is biologically represented in the panel

### macOS metadata files inside the proteome directory

They are ignored automatically, but you can clean them if you want:

```bash
find /path/to/proteomes -type f \( -name '._*' -o -name '.DS_Store' \) -delete
```

### codeml or HyPhy warnings

Warnings are tolerated when required result files are present. Hard failure occurs when required outputs are missing or unusable.

### `--resume` refuses my new analysis flags

This is intentional. Resume is designed to continue an existing run, not mutate its workflow definition midway through. Start a new `--outdir` if you want different analysis settings.

## Developer Notes

### Editable install

```bash
pip install -e .
```

### Build distributions

```bash
python -m pip install --upgrade build twine
python -m build --sdist --wheel
twine check dist/*
```

### Publish

```bash
twine upload dist/*
```

### Release checklist

1. Update version in `pyproject.toml` and `babappasnake/__init__.py`.
2. Run the test suite.
3. Build source and wheel distributions.
4. Check distributions.
5. Publish to PyPI.
6. Tag the release in Git.

## License

MIT
