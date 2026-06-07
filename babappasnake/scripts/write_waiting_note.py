#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path


def build_waiting_note(orthogroup: str, headers: str, expected_cds: str, outdir: str) -> str:
    expected_cds_path = Path(expected_cds).resolve()
    outdir_path = Path(outdir).resolve()
    return f"""babappasnake checkpoint reached.

Orthogroup proteins have been written to:
{orthogroup}

Protein headers for the orthogroup are listed in:
{headers}

When the corresponding CDS FASTA is ready, place it at:
{expected_cds_path}

Then resume from this checkpoint with:
babappasnake --resume --outdir {outdir_path}

Alternatively, let babappasnake stage the CDS file for you:
babappasnake --resume --outdir {outdir_path} --cds /path/to/orthogroup_cds.fasta

Important:
- CDS FASTA headers do not need to match protein headers.
- babappasnake will map CDS to orthogroup proteins by translated sequence similarity.
- Resume reloads the saved config.yaml, clears stale Snakemake locks, and continues from the CDS mapping stage instead of restarting orthogroup definition.
"""


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--orthogroup", required=True)
    p.add_argument("--headers", required=True)
    p.add_argument("--outfile", required=True)
    p.add_argument("--expected-cds", required=True)
    p.add_argument("--outdir", required=True)
    a = p.parse_args()
    text = build_waiting_note(a.orthogroup, a.headers, a.expected_cds, a.outdir)
    Path(a.outfile).write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
