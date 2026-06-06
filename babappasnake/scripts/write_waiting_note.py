#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--orthogroup", required=True)
    p.add_argument("--headers", required=True)
    p.add_argument("--outfile", required=True)
    p.add_argument("--expected-cds", required=True)
    a = p.parse_args()
    text = f"""babappasnake checkpoint reached.

One-to-one orthogroup proteins have been written to:
{a.orthogroup}

Protein headers for the orthogroup are listed in:
{a.headers}

Now download the CDS sequences corresponding to this orthogroup and save them as:
{a.expected_cds}

Important:
- CDS FASTA headers do not need to match protein headers.
- babappasnake will map CDS to orthogroup proteins by translated sequence similarity.
- After placing the CDS file, rerun the same babappasnake command.
"""
    Path(a.outfile).write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
