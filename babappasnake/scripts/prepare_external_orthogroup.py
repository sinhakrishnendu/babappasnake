#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from Bio import SeqIO


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--orthogroup-proteins", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--query-fasta", default="")
    return p.parse_args()


def first_query_id(query_fasta: str) -> str:
    if not query_fasta:
        return ""
    path = Path(query_fasta)
    if not path.exists() or path.stat().st_size == 0:
        return ""
    records = list(SeqIO.parse(str(path), "fasta"))
    return records[0].id if records else ""


def prepare_external_orthogroup(
    orthogroup_proteins: Path,
    outdir: Path,
    query_fasta: str = "",
) -> None:
    records = list(SeqIO.parse(str(orthogroup_proteins), "fasta"))
    if len(records) < 2:
        raise RuntimeError(
            "External orthogroup FASTA must contain at least two protein sequences "
            "(the query and at least one partner)."
        )

    ids = [record.id for record in records]
    duplicates = sorted({record_id for record_id in ids if ids.count(record_id) > 1})
    if duplicates:
        raise RuntimeError(
            "External orthogroup FASTA contains duplicate sequence IDs: "
            + ", ".join(duplicates)
        )

    query_id = first_query_id(query_fasta)
    if query_id not in ids:
        query_id = ids[0]

    outdir.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, outdir / "orthogroup_proteins.fasta", "fasta")
    (outdir / "orthogroup_headers.txt").write_text("\n".join(ids) + "\n", encoding="utf-8")

    rows = []
    partner_index = 0
    for record_id in ids:
        if record_id == query_id:
            continue
        partner_index += 1
        rows.append(
            {
                "species": f"external_{partner_index}",
                "query": query_id,
                "orthogroup": "external_curated",
                "orthology_mode": "external",
                "candidate_members": record_id,
                "n_candidate_members": "1",
                "selected_members": record_id,
                "n_selected_members": "1",
                "membership_class": "external_curated",
            }
        )

    with open(outdir / "orthogroup_summary.tsv", "w", encoding="utf-8", newline="") as fh:
        fieldnames = [
            "species",
            "query",
            "orthogroup",
            "orthology_mode",
            "candidate_members",
            "n_candidate_members",
            "selected_members",
            "n_selected_members",
            "membership_class",
        ]
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    metadata = {
        "source": "external",
        "query": query_id,
        "orthogroup": "external_curated",
        "retained_partner_count": len(rows),
        "retained_sequence_count_including_query": len(records),
        "summary": str((outdir / "orthogroup_summary.tsv").resolve()),
        "input_orthogroup_proteins": str(orthogroup_proteins.resolve()),
    }
    (outdir / "orthogroup_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    args = parse_args()
    prepare_external_orthogroup(
        orthogroup_proteins=Path(args.orthogroup_proteins),
        outdir=Path(args.outdir),
        query_fasta=str(args.query_fasta),
    )
