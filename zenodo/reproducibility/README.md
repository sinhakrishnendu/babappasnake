# Reproducibility Materials

This directory contains the validated example dataset and outputs used to confirm the `babappasnake` `v0.2.0` release.

## Included inputs

- `inputs/proteomes/`
  - cleaned proteome FASTA files used as `--input`
- `inputs/query_protein.fasta`
  - query protein FASTA used as `--query`
- `inputs/cds.fasta`
  - CDS FASTA used to seed the orthogroup CDS checkpoint with `--cds-input`

## Included outputs

- `results/results_pkg_live/`
  - curated results from the successful packaged run

## Replay command

See `run_metadata/validated_run_command.txt`.

The recommended archive-relative replay form is:

```bash
babappasnake \
  --input reproducibility/inputs/proteomes \
  --query reproducibility/inputs/query_protein.fasta \
  --outgroup "Culex quinquefasciatus" \
  --clipkit yes \
  --iqtree protein \
  --cds-input reproducibility/inputs/cds.fasta \
  --output reproduced_results \
  --threads 4 \
  --use-conda
```

This replay command assumes:

- the archive contents are the current working directory
- required external tools are available through the packaged Snakemake conda environments
