# BABAPPASnake
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19048331.svg)](https://doi.org/10.5281/zenodo.19048331)

`babappasnake` is a command-line workflow for episodic positive selection analysis on one orthogroup at a time. It is built for practical comparative genomics: reproducible run directories, guided stepwise execution, resumable restarts after interruption, and robustness summaries across multiple alignment and trimming pathways.

This README is the user manual for the current command-line interface (CLI) and workflow behavior.

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

1. Define an orthogroup from a protein query and proteome panel, or stage an externally curated orthogroup FASTA.
2. Map user-provided coding sequence (CDS) records to the selected orthogroup proteins with coding-quality filters.
3. Build protein and codon alignments across one or more alignment methods.
4. Run both `raw` and `clipkit` pathway variants for robustness.
5. Infer pathway-specific trees with IQ-TREE, or stage a user-supplied Newick tree for all pathways.
6. Optionally root trees with an outgroup text query.
7. Optionally screen for recombination with HyPhy genetic algorithm for recombination detection (GARD).
8. Run HyPhy adaptive branch-site random effects likelihood (aBSREL) and mixed effects model of evolution (MEME).
9. Select branch-site foregrounds dynamically from aBSREL output.
10. Run branch-site codeml.
11. Optionally run pathway-level codeml ancestral sequence reconstruction (ASR).
12. Optionally extract ancestor and descendant sequences plus branch substitutions.
13. Write pathway summaries, robustness reports, and run provenance.

Design assumptions:

- The workflow operates on one orthogroup at a time.
- CDS can be unavailable at the start. The pipeline supports a protein-first checkpointed workflow.
- Tree inference is optional. If no user tree is supplied, IQ-TREE builds pathway-specific trees. If a user tree is supplied, IQ-TREE tree inference is skipped and the supplied tree is reused downstream.
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

Always needed for a complete end-to-end run after orthogroup proteins are available:

- `hyphy`
- `codeml` from `paml`
- `clipkit`

Needed only when the workflow is asked to infer trees internally:

- `iqtree` or `iqtree2` or `iqtree3`

Needed only when using the default OrthoFinder-assisted orthogroup mode:

- `blastp`
- `makeblastdb`
- `orthofinder` for `--orthogroup-method orthofinder`

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
- The waiting note names the exact CDS entry point: `OUTDIR/user_supplied/orthogroup_cds.fasta`.
- You can then place CDS at that path and resume with `babappasnake --resume --outdir OUTDIR`, or stage a file with `babappasnake --resume --outdir OUTDIR --cds /path/to/orthogroup_cds.fasta`.

### `--orthogroup-proteins`

Optional externally curated orthogroup protein FASTA.

Behavior:

- When provided, BABAPPASNAKE skips OrthoFinder and starts from this protein FASTA.
- This is recommended when orthology has already been inferred, manually curated, or benchmarked outside the workflow.
- The FASTA must contain the query protein plus at least one partner sequence.
- If `--query` is also provided, its first sequence ID is used as the query ID in orthogroup metadata; otherwise the first record in `--orthogroup-proteins` is treated as the query.

### `--outgroup`

Optional text query used to root pathway trees by case-insensitive substring matching against tip labels.

Behavior:

- If omitted, rooting is skipped and the unrooted tree is propagated downstream.
- If provided but unmatched or too broad, rooting falls back safely to the unrooted tree instead of crashing the run.

### `--tree-mode` and `--tree`

Tree handling is controlled by `--tree-mode`.

Behavior:

- `--tree-mode iqtree` is the default and runs IQ-TREE for each selected method by trim pathway.
- `--tree-mode user --tree /path/to/tree.nwk` stages a user-supplied Newick tree and reuses it for every selected pathway.
- If a tree path is supplied with `--tree`, BABAPPASnake automatically treats the run as user-tree mode.
- In guided interactive mode, the workflow asks whether you want to supply a tree after CDS handling. Enter `no` or leave the tree prompt empty to let IQ-TREE build the trees.
- The user tree must use tip labels compatible with the orthogroup CDS/alignment labels used downstream.

## Orthogroup Discovery Modes

### External curated orthogroup: `--orthogroup-proteins`

Orthology inference can be handled outside BABAPPASNAKE.
Use `--orthogroup-proteins` to provide a curated protein FASTA from OrthoFinder, SonicParanoid, OMA, manual curation, or any other orthology workflow.
The pipeline then performs CDS mapping, alignment, tree inference, HyPhy/Phylogenetic Analysis by Maximum Likelihood (PAML) analyses, robustness summaries, and provenance reporting without running orthology inference internally.

### Default: `--orthogroup-method orthofinder`

OrthoFinder is the default built-in orthogroup helper when an external orthogroup FASTA is not supplied.

Behavior:

1. Run OrthoFinder.
2. BLAST the query against OrthoFinder orthogroup-member proteins.
3. Rank orthogroups by query support.
4. Extract members from the best supported orthogroup using `--orthology-mode`.

### OrthoFinder query mapping details

The workflow does not assume that the original query ID is directly present inside an OrthoFinder orthogroup. Instead it does a query-to-members mapping step:

1. Parse `Orthogroups.tsv`.
2. Build a combined FASTA of orthogroup-member proteins.
3. BLAST the query against that combined FASTA.
4. Filter by query coverage threshold.
5. Rank orthogroups by cross-species support and bitscore.
6. Extract the top supported orthogroup and apply the configured orthology/paralogy mode.

### Orthology and paralogy modes

- `--orthology-mode representative` is the default. Single-copy species are retained directly; multi-copy species contribute the best query-supported copy.
- `--orthology-mode strict` retains only species with exactly one member in the selected orthogroup.
- `--orthology-mode paralog` retains all copies from each species in the selected orthogroup.

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

- Orthogroup method defaults to `orthofinder`.
- Orthology mode defaults to `representative`.
- Alignment methods default to `4` which means all three aligners.
- Trimming is forced internally to robustness mode: `raw + clipkit`.
- Recombination defaults to `none`.
- ASR defaults to `yes`.
- CDS is requested only after orthogroup definition.
- Tree source is requested after CDS handling. Enter `no` if you want IQ-TREE to build pathway trees.
- Outgroup is requested after CDS handling and remains optional.

### Non-interactive batch mode

Use this for scripted or high-performance-computing-style runs:

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

### Example: representative orthology mode with ASR disabled

```bash
babappasnake \
  --prot /path/to/proteomes \
  --query /path/to/query.fasta \
  --cds /path/to/orthogroup_cds.fasta \
  --orthology-mode representative \
  --run-asr no \
  --outdir run_representative_no_asr \
  --threads 12 \
  --interactive no \
  --guided no
```

### Example: retain all paralog copies

```bash
babappasnake \
  --prot /path/to/proteomes \
  --query /path/to/query.fasta \
  --cds /path/to/orthogroup_cds.fasta \
  --orthology-mode paralog \
  --outdir run_paralog_mode \
  --threads 12 \
  --interactive no \
  --guided no
```

### Example: start from externally curated orthogroup proteins

```bash
babappasnake \
  --orthogroup-proteins /path/to/curated_orthogroup_proteins.fasta \
  --cds /path/to/orthogroup_cds.fasta \
  --outdir run_external_orthogroup \
  --threads 12 \
  --interactive no \
  --guided no
```

### Example: use a supplied phylogeny instead of IQ-TREE

```bash
babappasnake \
  --orthogroup-proteins /path/to/curated_orthogroup_proteins.fasta \
  --cds /path/to/orthogroup_cds.fasta \
  --tree-mode user \
  --tree /path/to/curated_species_tree.nwk \
  --outdir run_user_tree \
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

If the CDS file has already been placed at `run01/user_supplied/orthogroup_cds.fasta`, the shorter form is enough:

```bash
babappasnake --resume --outdir run01
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
- runs `snakemake --unlock` for that run directory before continuing
- calls Snakemake with `--rerun-incomplete`
- continues from missing or incomplete outputs instead of restarting from the initial entry point
- preserves completed outputs from earlier stages, including orthogroup definition after a CDS checkpoint

Guided resume behavior:

- saved and detected completed steps are skipped automatically
- the session restarts near the interruption point
- guided mode prints the exact `babappasnake --resume --outdir OUTDIR` command at the start and whenever the user stops at a step
- if CDS is not ready, guided mode names `OUTDIR/user_supplied/orthogroup_cds.fasta`, writes `WAITING_FOR_CDS.txt`, and does not ask downstream tree or outgroup questions until resume

Non-guided resume behavior:

- Snakemake reruns missing or incomplete work only
- if the run stopped at the CDS checkpoint, placing the CDS FASTA at `OUTDIR/user_supplied/orthogroup_cds.fasta` and running `babappasnake --resume --outdir OUTDIR` continues at `map_cds`

Files used for resume:

- `OUTDIR/config.yaml`
- `OUTDIR/.babappasnake/resume_state.json`

Allowed resume-time overrides:

- `--cds PATH`
- `--outgroup TEXT`
- `--tree PATH`
- `--tree-mode {iqtree,user}`
- `--threads INT`
- `--guided {yes,no}`
- `--snake-args "..." `

Analysis settings that define workflow structure are intentionally not changeable on resume. For example, do not expect `--resume` to accept a new orthogroup method, new alignment mode, new recombination mode, or new ASR setting for an existing run directory.

## Workflow Stages

1. `define_orthogroup`
   Build the selected orthogroup with OrthoFinder or stage externally curated orthogroup proteins.
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
7. `tree_all_pathways`
   Infer pathway trees with IQ-TREE or stage a user-supplied tree, depending on `--tree-mode`.
8. `root_iqtree_outgroup_all_pathways`
   Root trees when possible, otherwise propagate unrooted trees.
9. `hyphy_exploratory_all_pathways`
   Run aBSREL and MEME.
10. `parse_foregrounds_all_pathways`
    Select significant foreground branches from aBSREL.
11. `prepare_foreground_trees_all_pathways`
    Build branch-labeled trees for branch-site codeml.
12. `branchsite_batch_all_pathways`
    Run branch-site codeml and Benjamini-Hochberg (BH)-correct results.
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
- `orthogroup/orthogroup_summary.tsv`
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
- `--orthogroup-proteins PATH`
- `--outdir PATH`
- `--interactive {yes,no}`
- `--guided {yes,no}`
- `--resume`
- `--snake-args "..."`

### Orthogroup and alignment options

- `--orthogroup-method {orthofinder}`
- `--orthology-mode {strict,representative,paralog}`
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
- `--tree-mode {iqtree,user}`
- `--tree PATH`
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

1. Read `OUTDIR/orthogroup/WAITING_FOR_CDS.txt`.
2. Place the CDS FASTA at the exact path shown there, normally `OUTDIR/user_supplied/orthogroup_cds.fasta`.
3. Run `babappasnake --resume --outdir OUTDIR`.

You can also let the CLI stage the file:

```bash
babappasnake --resume --outdir OUTDIR --cds /path/to/orthogroup_cds.fasta
```

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

The workflow reloads saved settings, runs `snakemake --unlock`, and continues missing or incomplete work with `--rerun-incomplete`.

### OrthoFinder finds no usable orthogroup

The pipeline stops explicitly when the selected orthogroup and orthology mode retain no partner sequences.

Check:

- query quality
- proteome quality
- taxon sampling
- whether the query is biologically represented in the panel
- whether a curated external orthogroup should be supplied with `--orthogroup-proteins`

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
