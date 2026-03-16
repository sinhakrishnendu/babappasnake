#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path

from common import has_internal_stop_codon, read_fasta_records, write_fasta_records

STOP_CODONS = {"TAA", "TAG", "TGA"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint-dir", required=True)
    parser.add_argument("--validated-fasta-out", required=True)
    parser.add_argument("--member-check-out", required=True)
    parser.add_argument("--report-out", required=True)
    args = parser.parse_args()

    checkpoint_dir = Path(args.checkpoint_dir)
    metadata = json.loads((checkpoint_dir / "checkpoint_metadata.json").read_text())
    expected_cds = Path(metadata["expected_cds_fasta"])
    if not expected_cds.exists():
        raise FileNotFoundError(
            "CDS checkpoint reached. Provide the requested CDS FASTA and rerun Snakemake.\n"
            f"Expected file: {expected_cds}"
        )

    required_ids = metadata["required_member_ids"]
    provided_records = read_fasta_records(expected_cds)
    provided_by_id = {record["id"]: record for record in provided_records}
    missing = [record_id for record_id in required_ids if record_id not in provided_by_id]
    extra = [record["id"] for record in provided_records if record["id"] not in set(required_ids)]
    if missing or extra:
        raise ValueError(
            "The supplied CDS FASTA does not match the selected orthogroup members.\n"
            f"Missing IDs: {missing or 'none'}\n"
            f"Unexpected IDs: {extra or 'none'}"
        )

    validated_records = []
    member_rows = []
    trimmed_terminal_stops = 0
    for record_id in required_ids:
        sequence = provided_by_id[record_id]["seq"].upper().replace("U", "T")
        if any(base not in {"A", "C", "G", "T", "N"} for base in sequence):
            raise ValueError(f"Invalid nucleotide characters found in {record_id}")
        if len(sequence) >= 3 and sequence[-3:] in STOP_CODONS:
            sequence = sequence[:-3]
            trimmed_terminal_stops += 1
        if len(sequence) % 3 != 0:
            raise ValueError(f"Sequence length is not a multiple of three for {record_id}")
        if has_internal_stop_codon(sequence):
            raise ValueError(f"Internal stop codon detected in {record_id}")
        validated_records.append({"id": record_id, "seq": sequence})
        member_rows.append(
            {
                "member_id": record_id,
                "length_nt": len(sequence),
                "length_aa_estimate": len(sequence) // 3,
                "status": "validated",
            }
        )

    write_fasta_records(validated_records, args.validated_fasta_out)
    with open(args.member_check_out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(member_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(member_rows)

    Path(args.report_out).write_text(
        "\n".join(
            [
                f"Validated CDS FASTA: {args.validated_fasta_out}",
                f"Orthogroup members checked: {len(validated_records)}",
                f"Terminal stop codons trimmed: {trimmed_terminal_stops}",
                "Status: success",
            ]
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
