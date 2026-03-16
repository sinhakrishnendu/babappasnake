#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path

from common import has_internal_stop_codon, read_fasta_records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment-fasta", required=True)
    parser.add_argument("--members-tsv", required=True)
    parser.add_argument("--json-out", required=True)
    parser.add_argument("--report-out", required=True)
    args = parser.parse_args()

    alignment_records = read_fasta_records(args.alignment_fasta)
    if not alignment_records:
        raise ValueError("Aligned CDS FASTA is empty.")

    with open(args.members_tsv) as handle:
        expected_members = [row["member_id"] for row in csv.DictReader(handle, delimiter="\t")]

    observed_members = [record["id"] for record in alignment_records]
    if observed_members != expected_members:
        raise ValueError(
            "Aligned CDS sequence order does not match the selected orthogroup members.\n"
            f"Expected: {expected_members}\nObserved: {observed_members}"
        )

    alignment_length = len(alignment_records[0]["seq"])
    if alignment_length % 3 != 0:
        raise ValueError("Aligned CDS length is not a multiple of three.")

    sequence_summaries = []
    for record in alignment_records:
        if len(record["seq"]) != alignment_length:
            raise ValueError(f"Sequence {record['id']} has inconsistent alignment length.")
        ungapped = record["seq"].replace("-", "")
        if len(ungapped) % 3 != 0:
            raise ValueError(f"Ungapped CDS length is not a multiple of three for {record['id']}")
        if has_internal_stop_codon(record["seq"]):
            raise ValueError(f"Internal stop codon found after alignment in {record['id']}")
        gap_fraction = record["seq"].count("-") / alignment_length
        sequence_summaries.append(
            {
                "member_id": record["id"],
                "ungapped_length_nt": len(ungapped),
                "gap_fraction": round(gap_fraction, 4),
            }
        )

    validation = {
        "alignment_fasta": args.alignment_fasta,
        "sequence_count": len(alignment_records),
        "alignment_length_nt": alignment_length,
        "alignment_length_codons": alignment_length // 3,
        "sequence_summaries": sequence_summaries,
        "status": "ok",
    }
    Path(args.json_out).write_text(json.dumps(validation, indent=2))
    Path(args.report_out).write_text(
        "\n".join(
            [
                f"Alignment FASTA: {args.alignment_fasta}",
                f"Sequence count: {len(alignment_records)}",
                f"Aligned nucleotide length: {alignment_length}",
                "Validation: passed",
            ]
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
