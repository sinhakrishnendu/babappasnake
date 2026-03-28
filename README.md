# babappasnake
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19048331.svg)](https://doi.org/10.5281/zenodo.19048331)
`babappasnake` is a command-line workflow for episodic positive selection analysis on a single orthogroup.
It is designed for practical use: one command to launch, automatic checkpointing, resumable execution, and clear summary outputs.
In terminal mode, it can run as an interactive guided engine instead of a black-box single-shot command.

## What the pipeline does

1. Runs reciprocal best-hit (RBH) ortholog discovery.
2. Builds an orthogroup from your query and proteomes.
3. Maps user CDS records to orthogroup proteins (after lowercase intron clipping and uppercase ORF window extraction).
4. Creates protein and codon alignments for selected engines (`babappalign`, `mafft`, `prank`).
5. For `mafft`/`prank`, protein alignments are converted to codon alignments with a built-in robust Python back-translator (pal2nal-like behavior, tolerant to mild mismatches).
6. Branches each selected alignment method into trim pathways:
   - `raw` (untrimmed alignment branch)
   - `clipkit` (ClipKIT-trimmed branch)
7. Runs terminal stop-codon cleanup on codon alignments in both branches.
8. Infers an ML tree with IQ-TREE (`-m MFP -B 1000 -redo`) per `(method, trim_state)` pathway.
9. Roots inferred trees using a user-supplied outgroup label query (case-insensitive header matching; optional).
10. Runs HyPhy aBSREL and MEME with user-selected branch scopes (default `Leaves`) for each pathway.
11. Selects foreground branches from aBSREL using dynamic thresholding.
12. Runs branch-site `codeml` only for selected branches (alt and null models) per pathway.
13. Runs codeml ancestral sequence reconstruction (ASR) per pathway.
14. Produces pathway-level summaries plus robustness reports across both alignment method and trimming decision.

## Installation (for end users)

`pip` installs Python packages, but external bioinformatics binaries must also be available on `PATH`.

### Recommended setup (conda + pip)

```bash
conda create -n babappasnake -c conda-forge -c bioconda \
  python=3.11 blast iqtree hyphy paml clipkit pip
conda activate babappasnake
pip install babappasnake
```

Optional (only if you select those pathways):

```bash
conda install -c conda-forge -c bioconda mafft prank
```

Notes:
- `babappalign` is installed automatically as a PyPI dependency of `babappasnake`.
- `mafft`, `prank`, `blast`, `iqtree`, `hyphy`, and `paml/codeml` are external binaries and still need system/conda installation.

### Quick verification

```bash
babappasnake --help
which blastp makeblastdb hyphy codeml clipkit
which babappalign mafft prank
```

Notes:
- IQ-TREE binary detection is flexible (`iqtree2`, `iqtree3`, or `iqtree`).
- On Apple Silicon, `iqtree3` is common and is accepted automatically.

## Input requirements

- `--prot`: directory containing proteome FASTA files.
- `--query`: protein FASTA containing the query sequence.
- `--cds` (optional at first run): CDS FASTA for the orthogroup.
- `--outgroup`: outgroup query string used to root the IQ-TREE output (e.g., `culex` matches headers containing `culex`).
- `--alignment-methods`: numeric alignment selector (`1=babappalign`, `2=mafft`, `3=prank`, `4=all three`).

CDS quality checks (when `--cds` is supplied):
- Lowercase intron characters are clipped from each CDS.
- For each CDS, the best uppercase `ATG ... STOP` ORF window is retained.
- CDS records that still fail ORF/start-stop/frame checks are excluded.
- Proteins without a qualifying CDS match are skipped from codon/tree downstream analyses.

## Quick start

### Guided interactive mode (default in terminal)

```bash
babappasnake
```

This mode prompts for pipeline settings, executes one rule at a time, asks `run/skip/stop` at every step, and prints per-step output previews.
It asks for CDS only after `rbh_orthogroup` finishes, then asks optional outgroup text for rooting.
It also prints explicit orthogroup membership in terminal: groups included and groups omitted at RBH stage.
If outgroup is left empty, rooting is safely skippable in guided mode and downstream uses unrooted IQ-TREE trees.
You can choose one method or all three methods from the numbered selector. Default is `4` (all three).
You can choose trimming strategy:
- `raw` only
- `clipkit` only
- `both` (robustness mode; runs both branches for each selected method)
Available cores are split across active pathways (`selected_methods x selected_trim_states`), and total cores are auto-raised only when below pathway count.

### Case A: you already have the CDS file

```bash
babappasnake \
  --prot /path/to/proteomes \
  --query /path/to/query.fasta \
  --cds /path/to/orthogroup_cds.fasta \
  --alignment-methods 4 \
  --trim-strategy both \
  --outgroup culex \
  --outdir run01 \
  --threads 12 \
  --interactive no \
  --guided no
```

### Case B: two-stage run (CDS provided later)

Run once:

```bash
babappasnake \
  --prot /path/to/proteomes \
  --query /path/to/query.fasta \
  --outgroup culex \
  --outdir run01 \
  --threads 12 \
  --interactive no \
  --guided no
```

The first run stops intentionally and writes:

- `run01/orthogroup/orthogroup_proteins.fasta`
- `run01/orthogroup/orthogroup_headers.txt`
- `run01/orthogroup/WAITING_FOR_CDS.txt`

Then place your CDS FASTA at:

- `run01/user_supplied/orthogroup_cds.fasta`

Re-run the same command. The workflow resumes automatically.

## Robustness pathway model

By design, the current pipeline always evaluates both trim states (`raw` and `clipkit`) for each selected alignment method.  
For `--alignment-methods 4`, one run evaluates:

- `babappalign_raw`
- `babappalign_clipkit`
- `mafft_raw`
- `mafft_clipkit`
- `prank_raw`
- `prank_clipkit`

Each pathway keeps isolated outputs under `<module>/<method>/<trim_state>/...` to avoid collisions.

## Dynamic significance logic

Foreground selection from aBSREL uses dynamic p-thresholding:

- start at `p <= 0.05`
- if no branch passes, increase by `0.01`
- stop as soon as at least one branch is found
- hard upper bound: `0.2`

Only those selected branches go to branch-site `codeml`.
Each selected branch runs two codeml fits (alternative + null), then BH-FDR correction is applied.

## Output guide

All outputs are written under `--outdir` (default: `babappasnake_run`).

Most important files:

- `summary/episodic_selection_summary.txt`: top-level summary alias (primary selected method).
- `summary/<method>/<trim_state>/episodic_selection_summary.txt`: pathway-specific final reports.
- `summary/robustness_matrix.tsv`: one row per `(method, trim_state)` pathway with alignment/tree/test counts and status.
- `summary/robustness_consensus.tsv`: replicated signals across pathways with reproducibility labels.
- `summary/robustness_narrative.txt`: human-readable robustness interpretation.
- `summary/comparative_reproducibility_summary.txt`: alignment-method sensitivity, trim sensitivity, and label counts.
- `summary/robustness_publication_table.tex`: publication-ready comparative signal table.
- `summary/run_provenance.json`: machine-readable run provenance (methods, trim states, parameters, key outputs).
- `hyphy/<method>/<trim_state>/foreground_threshold.json`: selected dynamic threshold and hit count.
- `hyphy/<method>/<trim_state>/significant_foregrounds.tsv`: selected aBSREL foreground branches.
- `branchsite/<method>/<trim_state>/branchsite_results.tsv`: codeml branch-site statistics and BH significance.
- `asr/<method>/<trim_state>/asr_done.json`: pathway ASR completion record.
- `tree/<method>/<trim_state>/orthogroup.treefile`: pathway inferred ML tree (unrooted IQ-TREE output).
- `tree/<method>/<trim_state>/orthogroup.rooted.treefile`: pathway rooted tree used downstream.

## CLI reference

```text
babappasnake --prot PROTEOMES_DIR --query QUERY_FASTA [options]
```

Options:

- `--cds PATH`: CDS FASTA (optional for initial run).
- `--outgroup TEXT`: outgroup query used for tree rooting (case-insensitive substring match against tip headers).
- `--outdir PATH`: output directory (default: `babappasnake_run`).
- `--alignment-methods {1,2,3,4}`: method selection (`1=babappalign`, `2=mafft`, `3=prank`, `4=all three`; default: `4`).
- `--trim-strategy {raw,clipkit,both}`: accepted for compatibility, but runtime is forced to `both` so both raw and clipkit summaries are always generated.
- `--coverage FLOAT`: RBH reciprocal coverage minimum (default: `0.70`).
- `--threads INT`: total Snakemake cores. Active pathways receive an equal split (`floor(threads / (selected_methods x selected_trim_states))`), and total cores are auto-raised only when below pathway count (default: detected CPU core count).
- `--iqtree-bootstrap INT`: UFBoot replicates for IQ-TREE (default: `1000`; typical options: `1000`, `5000`, `10000`).
- `--iqtree-bnni {yes,no}`: enable/disable IQ-TREE `-bnni` (default: `no`).
- `--iqtree-model TEXT`: IQ-TREE model string (default: `MFP`).
- `--absrel-branches TEXT`: HyPhy aBSREL branch selector (default: `Leaves`; common choices: `Leaves`, `Internal`, `All`).
- `--meme-branches TEXT`: HyPhy MEME branch selector (default: `Leaves`; common choices: `Leaves`, `Internal`, `All`).
- `--codeml-codonfreq INT`: codeml `CodonFreq` value for branch-site and ASR runs (default: `7`; e.g. `1`, `2`, `7`).
- `--absrel-p FLOAT`: compatibility parameter retained in config (default: `0.05`; dynamic mode is used for branch selection).
- `--absrel-dynamic-start FLOAT`: dynamic foreground start p-value (default: `0.05`).
- `--absrel-dynamic-step FLOAT`: dynamic foreground increment (default: `0.01`).
- `--absrel-dynamic-max FLOAT`: dynamic foreground max p-value (default: `0.2`).
- `--meme-p FLOAT`: MEME reporting threshold in summary (default: `0.05`).
- `--use-clipkit {yes,no}`: legacy compatibility argument; current robustness mode always runs both raw and clipkit branches.
- `--clipkit-mode-protein TEXT`: ClipKIT mode for protein trimming (default: `kpic-smart-gap`).
- `--clipkit-mode-codon TEXT`: ClipKIT mode for codon trimming (default: `kpic-smart-gap`).
- `--snake-args "..."`: extra raw arguments forwarded to Snakemake.
- `--interactive {yes,no}`: prompt for settings at launch (default: `yes`).
- `--guided {yes,no}`: execute one rule at a time with confirmation (default: `yes`).

Example with additional Snakemake flags:

```bash
babappasnake \
  --prot /path/to/proteomes \
  --query /path/to/query.fasta \
  --cds /path/to/orthogroup_cds.fasta \
  --alignment-methods 4 \
  --trim-strategy both \
  --outgroup culex \
  --outdir run01 \
  --threads 12 \
  --snake-args "--keep-going"
```

## Troubleshooting

### "Missing required external tools"

Install missing binaries and ensure they are on `PATH` in the same shell where you run `babappasnake`.

### Run stops with `WAITING_FOR_CDS.txt`

This is expected for two-stage usage.
Add `user_supplied/orthogroup_cds.fasta` and re-run.

### codeml returns non-zero or writes warnings

The workflow accepts codeml warning-heavy runs when valid output files are produced.
Hard failure is triggered only when required codeml result files are missing.

### Can I resume after interruption?

Yes. Re-run the same command; Snakemake resumes and re-runs incomplete jobs as needed.

## Developer/local source install

For local development:

```bash
pip install -e .
```

The package entry-point command is still `babappasnake`.

## Maintainer release checklist (GitHub + PyPI)

1. Update version in `pyproject.toml`.
2. Commit and push to GitHub.
3. Build distributions:

```bash
python -m pip install --upgrade build twine
python -m build
twine check dist/*
```

4. Publish to PyPI:

```bash
twine upload dist/*
```

5. Optionally create and push a matching git tag.

## License

MIT
