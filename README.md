# BABAPPASnake
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19048331.svg)](https://doi.org/10.5281/zenodo.19048331)

`babappasnake` is a reproducible command-line workflow for episodic positive selection analysis on one orthogroup at a time.
It is designed for practical comparative genomics: resumable runs, interactive stepwise control, and robustness summaries across alignment and trimming choices.

## Quick Start

### 1) Install tools and package

```bash
conda create -n babappasnake -c conda-forge -c bioconda \
  python=3.11 blast orthofinder iqtree hyphy paml clipkit mafft prank pip
conda activate babappasnake
pip install babappasnake
```

### 2) Run (non-interactive)

```bash
babappasnake \
  --prot /path/to/proteomes \
  --query /path/to/query.fasta \
  --cds /path/to/orthogroup_cds.fasta \
  --alignment-methods 4 \
  --outgroup culex \
  --outdir run01 \
  --threads 12 \
  --interactive no \
  --guided no
```

### 3) Run (interactive guided mode)

```bash
babappasnake
```

This prompts step by step, shows what each stage does, and supports `run/skip/stop` per stage.

## Manual Contents

1. [What The Pipeline Does](#what-the-pipeline-does)
2. [Installation And Environment](#installation-and-environment)
3. [Input Requirements](#input-requirements)
4. [Orthogroup Discovery Strategy](#orthogroup-discovery-strategy)
5. [Running Modes](#running-modes)
6. [CDS Mapping And QC](#cds-mapping-and-qc)
7. [Alignment, Trimming, Tree, And Selection Steps](#alignment-trimming-tree-and-selection-steps)
8. [ASR Extraction Of Selected Branches](#asr-extraction-of-selected-branches)
9. [Outputs And Directory Layout](#outputs-and-directory-layout)
10. [Resume, Rerun, And Reproducibility](#resume-rerun-and-reproducibility)
11. [CLI Reference](#cli-reference)
12. [Troubleshooting](#troubleshooting)
13. [Developer And Release Notes](#developer-and-release-notes)

## What The Pipeline Does

A complete run performs these stages:

1. Build orthogroup candidates from your query and proteomes.
2. Select the best orthogroup strategy by strict 1:1 ortholog support.
3. Map user CDS to selected proteins with coding-quality filtering.
4. Build protein alignments with selected MSA engines.
5. Build codon alignments (native for BABAPPAlign; robust back-translation for MAFFT/PRANK).
6. Run both trimming states (`raw` and `clipkit`) for robustness.
7. Infer trees with IQ-TREE for each `(method, trim_state)` pathway.
8. Optionally root trees with outgroup text query.
9. Optionally run HyPhy GARD recombination screening per pathway.
10. Run HyPhy aBSREL + MEME per pathway.
11. Select foreground branches dynamically from aBSREL.
12. Run branch-site codeml per selected foreground.
13. Run codeml ASR per pathway.
14. Extract ancestor/descendant sequences and substitutions for selected branches.
15. Write pathway summaries plus cross-pathway robustness reports.

## Installation And Environment

### Required external tools

- `blastp`
- `makeblastdb`
- `orthofinder`
- `iqtree` (or `iqtree2` or `iqtree3`)
- `hyphy`
- `codeml` (from `paml`)
- `clipkit`

### Optional but strongly recommended

- `mafft` (if using alignment method 2 or 4)
- `prank` (if using alignment method 3 or 4)

### Python package

`babappalign` is installed automatically as a Python dependency of `babappasnake`.

### Verify installation

```bash
babappasnake --help
which blastp makeblastdb orthofinder iqtree iqtree2 iqtree3 hyphy codeml clipkit
which babappalign mafft prank
```

## Input Requirements

### `--prot` proteomes directory

- Directory of protein FASTA files, one file per species.
- Supported extensions: `.fa`, `.faa`, `.fasta` (case-insensitive).
- Hidden/macOS metadata files are ignored (e.g. `.DS_Store`, `._*`, hidden files).
- Empty or malformed FASTA files are skipped with warnings.

### `--query` query FASTA

- Protein FASTA with exactly one query sequence.

### `--cds` CDS FASTA

- Optional on first run.
- If not provided initially, workflow stops after orthogroup definition and writes `WAITING_FOR_CDS.txt`.
- Add CDS file at `OUTDIR/user_supplied/orthogroup_cds.fasta` and rerun.

### `--outgroup` optional

- String used to root tree by case-insensitive substring match on tip headers.
- If omitted, unrooted trees are used downstream.

## Orthogroup Discovery Strategy

### Default strategy (`--orthogroup-method rbh`)

In default mode, orthogroup selection uses RBH only:

1. Run RBH stage.
2. Compute strict 1:1 ortholog count from RBH output.
3. If RBH yields zero strict 1:1 orthologs, stop with explicit error.

### Fallback comparison mode (`--orthogroup-method rbh_fallback`)

This preserves the previous behavior as an explicit opt-in:

1. Run RBH stage.
2. Run OrthoFinder stage.
3. Compute strict 1:1 ortholog count for each.
4. Select backend with larger strict 1:1 count.
5. If tied, keep RBH deterministically.
6. If both have zero strict 1:1 orthologs, stop with explicit error.

The selected backend and counts are printed explicitly.

### How OrthoFinder query mapping is done

OrthoFinder query mapping is BLAST-based, not query-ID membership based:

1. Parse OrthoFinder `Orthogroups.tsv` to load all groups and members.
2. Build a combined FASTA of orthogroup-member proteins with subject IDs encoded as:
   `<orthogroup>||<species>||<member>`
3. Run `blastp(query -> combined orthogroup members)`.
4. Filter by query coverage threshold.
5. Rank orthogroups by:
   - number of species with passing hits,
   - summed best bitscore per species,
   - top bitscore,
   - orthogroup ID (stable tie-break).
6. Select top-ranked orthogroup for downstream extraction.

### Strict 1:1 rule used downstream

A species contributes only if exactly one ortholog is retained for that species.
This guarantees no duplicate ortholog entries for the same species in the final selected orthogroup.

### Direct OrthoFinder mode (`--orthogroup-method orthofinder`)

- Runs OrthoFinder selection directly.
- Uses the same BLAST-based query-to-orthogroup mapping and strict 1:1 filtering.

## Running Modes

## Interactive guided mode (default)

```bash
babappasnake
```

Behavior:

- Prompts for required settings.
- Executes one rule at a time.
- Asks `run/skip/stop` for each stage.
- Shows step outputs and previews.
- Auto-skips already-completed steps safely.

Important guided-mode defaults:

- Orthogroup backend defaults to `rbh` in interactive mode.
- You can choose `rbh`, `orthofinder`, or `rbh_fallback`.
- Trimming is forced to robustness mode (`raw + clipkit`) for comparative summaries.
- CDS is asked only after orthogroup stage finishes.
- Outgroup prompt comes after CDS prompt and is optional.

## Non-interactive mode

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

Use this for scripted runs and HPC wrappers.

### Example: enable optional GARD screening

```bash
babappasnake \
  --prot /path/to/proteomes \
  --query /path/to/query.fasta \
  --cds /path/to/orthogroup_cds.fasta \
  --alignment-methods 4 \
  --recombination gard \
  --gard-mode Faster \
  --gard-rate-classes 3 \
  --outdir run01_gard \
  --threads 12 \
  --interactive no \
  --guided no
```

## Two-stage run (no CDS at start)

Stage 1:

```bash
babappasnake \
  --prot /path/to/proteomes \
  --query /path/to/query.fasta \
  --outdir run01 \
  --interactive no \
  --guided no
```

After stage 1 completes, add:

- `run01/user_supplied/orthogroup_cds.fasta`

Then rerun same command to resume.

## CDS Mapping And QC

During CDS mapping:

1. Lowercase intronic segments are clipped out.
2. Uppercase ORF window is retained.
3. Must start with uppercase start codon and end with uppercase stop codon.
4. Frame consistency checks are applied.
5. Failing CDS entries are excluded with warnings.

Outputs:

- `mapped_cds/mapped_orthogroup_cds.fasta`
- `mapped_cds/mapped_orthogroup_proteins.fasta`
- `mapped_cds/cds_protein_mapping.tsv`

## Alignment, Trimming, Tree, And Selection Steps

## Alignment methods

`--alignment-methods` options:

- `1`: `babappalign`
- `2`: `mafft`
- `3`: `prank`
- `4`: all three

## Trimming model

Pipeline is enforced to run both trim states for robustness:

- `raw`
- `clipkit`

So each selected method is expanded into two pathways.
Example for method `4`:

- `babappalign_raw`
- `babappalign_clipkit`
- `mafft_raw`
- `mafft_clipkit`
- `prank_raw`
- `prank_clipkit`

## Tree inference and rooting

- IQ-TREE runs per pathway.
- Outgroup rooting is optional and applied if query text is supplied and matched.
- If no outgroup is supplied/matched, downstream continues with unrooted tree.

## Optional recombination screening (HyPhy GARD)

- Controlled by `--recombination {none,gard,auto}` (default: `none`).
- `none`: skip recombination screening (legacy/default behavior).
- `gard`: run HyPhy GARD per active `(method, trim_state)` pathway.
- `auto`: currently an alias of `gard`.
- Additional controls:
  - `--gard-mode {Normal,Faster}` (default: `Faster`)
  - `--gard-rate-classes INT` (default: `3`)

Output locations:

- `recombination/<method>/<trim_state>/gard/gard.json`
- `recombination/<method>/<trim_state>/gard/gard_summary.json`
- `recombination/<method>/<trim_state>/gard/gard.stdout.txt`
- `recombination/<method>/<trim_state>/gard/gard.stderr.txt`

Interpretation note:

- In the current implementation, GARD is a screening/reporting module.
- Downstream branch-site/codeml analyses remain full-length by default unless explicit fragment-aware routing is added in a future release.

## HyPhy and branch-site selection

- aBSREL and MEME run per pathway.
- Foregrounds are selected by dynamic aBSREL threshold:
  - start `0.05`
  - increment `0.01`
  - cap `0.2`
- Selected foregrounds feed branch-site codeml.
- Branch-site outputs are BH-corrected.

## ASR Extraction Of Selected Branches

After branch-site selection, `extract_selected_branch_ancestors` does:

1. Map selected branches onto canonical tree edges `(parent_node -> child_node)`.
2. Recover ancestor and descendant CDS/AA sequences.
3. Compute codon and amino-acid substitutions per selected branch.
4. Annotate overlaps with MEME/BEB where available.
5. Write branch-level summary and provenance.

This stage is model-based reconstruction from codeml outputs.

## Outputs And Directory Layout

All outputs are inside `--outdir`.

High-value files:

- `orthogroup/orthogroup_proteins.fasta`
- `orthogroup/orthogroup_headers.txt`
- `orthogroup/rbh_summary.tsv`
- `mapped_cds/cds_protein_mapping.tsv`
- `alignments/<method>/<trim_state>/...`
- `tree/<method>/<trim_state>/orthogroup.treefile`
- `tree/<method>/<trim_state>/orthogroup.rooted.treefile`
- `hyphy/<method>/<trim_state>/absrel.json`
- `hyphy/<method>/<trim_state>/meme.json`
- `hyphy/<method>/<trim_state>/significant_foregrounds.tsv`
- `recombination/<method>/<trim_state>/gard/gard_summary.json`
- `branchsite/<method>/<trim_state>/branchsite_results.tsv`
- `asr/<method>/<trim_state>/asr_done.json`
- `asr/branch_to_nodes.tsv`
- `asr/branch_substitutions.tsv`
- `asr/selected_branch_asr_summary.tsv`
- `summary/<method>/<trim_state>/episodic_selection_summary.txt`
- `summary/robustness_matrix.tsv`
- `summary/robustness_consensus.tsv`
- `summary/robustness_narrative.txt`
- `summary/comparative_reproducibility_summary.txt`
- `summary/robustness_publication_table.tex`
- `summary/run_provenance.json`

Top-level compatibility aliases:

- `summary/episodic_selection_summary.txt`
- `asr/asr_done.json`

## Resume, Rerun, And Reproducibility

- Re-running the same command resumes existing work.
- In guided mode, completed steps are auto-detected and skipped.
- Final provenance is written to `summary/run_provenance.json`.
- ASR extraction provenance is written to `asr/asr_extraction_provenance.json`.

## CLI Reference

Basic form:

```text
babappasnake --prot PROTEOMES_DIR --query QUERY_FASTA [options]
```

Core options:

- `--cds PATH`
- `--orthogroup-method {rbh,orthofinder,rbh_fallback}`
- `--alignment-methods {1,2,3,4}`
- `--outgroup TEXT`
- `--outdir PATH`
- `--threads INT`
- `--interactive {yes,no}`
- `--guided {yes,no}`

Selection and model options:

- `--coverage FLOAT` (RBH/BLAST mapping coverage threshold)
- `--iqtree-bootstrap INT`
- `--iqtree-bnni {yes,no}`
- `--iqtree-model TEXT`
- `--absrel-branches TEXT`
- `--meme-branches TEXT`
- `--codeml-codonfreq INT`
- `--recombination {none,gard,auto}`
- `--gard-mode {Normal,Faster}`
- `--gard-rate-classes INT`
- `--absrel-p FLOAT`
- `--absrel-dynamic-start FLOAT`
- `--absrel-dynamic-step FLOAT`
- `--absrel-dynamic-max FLOAT`
- `--meme-p FLOAT`
- `--clipkit-mode-protein TEXT`
- `--clipkit-mode-codon TEXT`
- `--snake-args "..."`

Compatibility options:

- `--trim-strategy {raw,clipkit,both}` accepted, but runtime robustness mode is enforced to `both`.
- `--use-clipkit {yes,no}` retained for backward compatibility.

## Troubleshooting

## "Missing required external tools"

Install missing binaries in the active environment and verify with `which`.

## Run stops at `WAITING_FOR_CDS.txt`

This is expected when CDS is not yet supplied.
Add `user_supplied/orthogroup_cds.fasta` and rerun.

## OrthoFinder/RBH finds no usable orthogroup

Pipeline stops explicitly if both strategies yield zero strict 1:1 ortholog support.
Check proteome quality, query quality, and species composition.

## macOS metadata files in proteomes

You can clean sidecar files before running:

```bash
find /path/to/proteomes -type f \( -name '._*' -o -name '.DS_Store' \) -delete
```

## codeml warnings

Warnings are tolerated when required output files are present.
Hard failure occurs only when mandatory result files are missing.

## Developer And Release Notes

## Local editable install

```bash
pip install -e .
```

## Build package

```bash
python -m pip install --upgrade build twine
python -m build --sdist --wheel
twine check dist/*
```

## Publish to PyPI

```bash
twine upload dist/*
```

## Release checklist

1. Update version in `pyproject.toml`.
2. Run tests.
3. Build distributions.
4. Publish to PyPI.
5. Push Git tag and GitHub release.

## License

MIT
