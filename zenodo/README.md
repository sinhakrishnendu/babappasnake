# Zenodo Archive Bundle for babappasnake v0.2.0

This directory is a self-contained archival bundle prepared for upload to the Zenodo record targeted by:

- `10.5281/zenodo.19045558`

The bundle contains the exact `v0.2.0` release snapshot used after the codeml bug-fix validation pass, plus the validated real-data example run needed for later reuse and audit.

## Contents

- `source/`
  - exact `v0.2.0` source snapshot required to rebuild or inspect the workflow
- `distributions/`
  - built wheel and source distribution artifacts for `babappasnake==0.2.0`
- `reproducibility/inputs/`
  - clean proteome FASTAs plus the query protein and CDS input used in the validated run
- `reproducibility/results/`
  - curated copy of the validated `results_pkg_live` run outputs
- `reproducibility/run_metadata/`
  - runtime configuration, software versions, path remapping notes, and a replay command
- `SHA256SUMS.txt`
  - file checksums for the archive payload
- `CONTENTS.tsv`
  - top-level manifest of the archived materials

## Notes

- Local machine-only files were excluded from this archive, including local conda/keyring files and the path-specific `config/package_live_override.yaml`.
- Stale internal-node codeml artifacts from the earlier pre-fix run were excluded from the archived results. The archived codeml results only represent the final leaf-only policy.
- Some raw workflow outputs preserve absolute paths from the original validation machine. See `reproducibility/run_metadata/path_remapping.txt` for the archive-relative mapping.

## Primary final results

- `reproducibility/results/results_pkg_live/reports/executive_summary.txt`
- `reproducibility/results/results_pkg_live/codeml/codeml_branchsite_summary.tsv`
- `reproducibility/results/results_pkg_live/codeml/codeml_branchsite_summary.json`
