# Changelog

## v0.2.0

- Restricted codeml foreground selection to terminal aBSREL branches only; internal `Node*` branches are excluded from both significant-branch and fallback selection paths.
- Made codeml branch-site execution robust to the known non-fatal `codeml` exit-status `1` behavior by accepting runs that still produce parseable null and alternative lnL outputs.
- Preserved explicit null-model and alternative branch-site outputs for every tested branch and reported `lnL_alt`, `lnL_null`, LRT, raw p-values, and BH-adjusted significance in the executive summary.
- Validated the packaged workflow end-to-end on the local real dataset and generated a curated Zenodo reproducibility bundle under `zenodo/`.
