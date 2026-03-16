# BABAPPAsnake

`babappasnake` is a packaged Snakemake workflow for exploratory orthology-based molecular evolution analysis. It starts from a directory of proteomes plus a query protein, finds RBH-based ortholog candidates, selects a one-to-one orthogroup with mandatory outgroup retention, pauses at a resumable CDS checkpoint, and then runs BABAPPAlign, optional ClipKIT trimming, IQ-TREE 3, HyPhy, selective codeml branch-site tests, BH correction, and executive reporting.

## Release

Current release prepared in this repository: `v0.2.0`

`v0.2.0` includes the validated codeml bug-fix release used for the packaged end-to-end real-data run:

- codeml foregrounds are restricted to terminal aBSREL branches only
- codeml runs both branch-site alternative and null models for every tested branch
- LRT, raw p-values, and BH-adjusted significance are reported explicitly
- non-fatal `codeml` exit-status `1` behavior is tolerated when parseable outputs are present

## Installation

For a local editable install from this repository:

```bash
pip install -e .
```

For a regular install after publishing:

```bash
pip install babappasnake
```

Confirm the CLI is available:

```bash
babappasnake --help
```

## Direct CLI Usage

The direct run form is:

```bash
babappasnake \
  --input "/path/to/proteomes" \
  --query "/path/to/query.fasta" \
  --outgroup "Culex quinquefasciatus" \
  --clipkit yes \
  --iqtree protein
```

Alternative direct run form with trimming disabled and IQ-TREE on CDS:

```bash
babappasnake \
  --input "/path/to/proteomes" \
  --query "/path/to/query.fasta" \
  --outgroup "Culex quinquefasciatus" \
  --clipkit no \
  --iqtree cds
```

Required user-facing arguments:

- `--input`: path to the proteome directory
- `--query`: query protein FASTA
- `--outgroup`: exact outgroup taxon name
- `--clipkit`: `yes` or `no`
- `--iqtree`: `protein` or `cds`

Common optional arguments:

- `--config /path/to/config.yaml`: use a base YAML config and let the CLI override the main run settings
- `--threads 8`: Snakemake core count
- `--output results_run1`: results directory
- `--coverage-thresholds 50,60,70,80,90,95`
- `--cds-input /path/to/cds.fasta`: seed the checkpoint with an existing CDS file
- `--use-conda`
- `--dry-run`

Write the bundled default config template without running the workflow:

```bash
babappasnake --init-config babappasnake.config.yaml
```

Validate that the installed package can locate its bundled workflow:

```bash
babappasnake --validate-installation
```

Print the executive summary from a finished run:

```bash
babappasnake --summarize /path/to/results
```

## Workflow Routing Logic

BABAPPAlign always produces both of these files after the CDS checkpoint:

- `results/alignment/babappalign/protein/aligned_proteins.fasta`
- `results/alignment/babappalign/cds/aligned_cds.fasta`

If `--clipkit yes`, the workflow runs ClipKIT in `kpic-smartgap` mode on the protein alignment and keeps:

- `results/alignment/clipkit/protein/trimmed_proteins.fasta`
- `results/alignment/clipkit/cds/trimmed_cds.fasta`
- `results/alignment/clipkit/cds/projected_trimmed_cds.fasta`

`trimmed_cds.fasta` is the default biologically authoritative codon-safe CDS alignment. It is produced directly by ClipKIT in codon mode and is used downstream everywhere whenever `--clipkit yes`.

`projected_trimmed_cds.fasta` is preserved as a protein-guided QC/audit artifact so the retained protein mask can still be compared against the direct codon-mode CDS trim.

If `--clipkit no`, trimming is skipped and downstream steps use the untrimmed BABAPPAlign alignments.

IQ-TREE input selection is explicit:

- `--iqtree protein`
  - with `--clipkit yes`: IQ-TREE uses `trimmed_proteins.fasta`
  - with `--clipkit no`: IQ-TREE uses `aligned_proteins.fasta`
- `--iqtree cds`
  - with `--clipkit yes`: IQ-TREE uses codon-safe `trimmed_cds.fasta`
  - with `--clipkit no`: IQ-TREE uses `aligned_cds.fasta`

HyPhy and codeml always use the codon-compatible CDS alignment:

- `trimmed_cds.fasta` when `--clipkit yes`
- `aligned_cds.fasta` when `--clipkit no`

The tree used by HyPhy and codeml is whichever rooted IQ-TREE tree was produced from the user-selected `--iqtree` mode.

For codeml branch-site follow-up, only terminal branches from the aBSREL exploratory run are considered. Internal `Node*` branches are excluded from codeml target selection, including the fallback lowest-p branch logic.

For every selected codeml branch, the workflow runs both:

- the branch-site alternative model
- the branch-site null model

The final codeml report is derived from the likelihood-ratio test `2 * (lnL_alt - lnL_null)`, followed by BH correction across tested branches.

## Checkpoint And Resume

The workflow intentionally pauses after orthogroup selection. It writes a checkpoint request directory containing the required member IDs and a README describing the expected CDS input.

When you reach the checkpoint, place the CDS FASTA at:

- `results/cds_input_checkpoint/request/user_supplied_cds.fasta`

Then rerun the exact same `babappasnake ...` command. Snakemake resumes from the checkpoint.

If you already have the CDS FASTA before the run starts, pass it with `--cds-input` so the checkpoint is pre-seeded automatically.

## Output Structure

```text
results/
  rbh/
  threshold_comparison/
  selected_orthogroup/
  cds_input_checkpoint/
    request/
    validated/
  alignment/
    babappalign/
      protein/
      cds/
    clipkit/
      protein/
      cds/
  iqtree/
  hyphy/
    absrel/
    busted/
    meme/
  codeml/
  reports/
  logs/
```

Key outputs:

- `results/selected_orthogroup/selection_report.txt`
- `results/alignment/alignment_validation.json`
- `results/iqtree/rooted_labeled.treefile`
- `results/hyphy/absrel/absrel_branch_summary.tsv`
- `results/hyphy/busted/busted_summary.json`
- `results/hyphy/meme/meme_sites.tsv`
- `results/codeml/codeml_branchsite_summary.tsv`
- `results/reports/executive_summary.txt`

## Configuration

The packaged default config lives at `config/config.yaml` in the source tree and is bundled into the installed package. The CLI writes a merged runtime copy under:

- `<results_dir>/.babappasnake/runtime_config.yaml`

You can still run the workflow directly with Snakemake if needed:

```bash
python -m snakemake \
  --snakefile Snakefile \
  --configfile config/config.yaml \
  --cores 8 \
  --rerun-incomplete
```

## Development Notes

Build a distributable package locally:

```bash
python -m build
```

Run the tests:

```bash
python -m pytest -q
```

The package build refreshes the bundled workflow assets automatically so the installed `babappasnake` command can launch Snakemake without needing the source tree.

## Zenodo Bundle

A curated reproducibility bundle for archival upload is written under `zenodo/` in this repository. It contains:

- the `v0.2.0` source snapshot
- built wheel and source distribution artifacts
- the validated example input dataset used for the packaged run
- the validated `results_pkg_live` workflow outputs
- release metadata, checksums, and run-manifest files
