#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

from common import read_fasta_records, write_fasta_records


def require_alignment(records, label):
    if not records:
        raise ValueError(f"{label} is empty.")
    expected_length = len(records[0]["seq"])
    ids = [record["id"] for record in records]
    if len(set(ids)) != len(ids):
        raise ValueError(f"{label} contains duplicate sequence identifiers.")
    for record in records:
        if len(record["seq"]) != expected_length:
            raise ValueError(f"{label} contains inconsistent alignment lengths.")
    return expected_length


def require_same_ids(reference_records, other_records, reference_label, other_label):
    reference_ids = [record["id"] for record in reference_records]
    other_ids = [record["id"] for record in other_records]
    if reference_ids != other_ids:
        raise ValueError(
            f"{other_label} sequence IDs do not match {reference_label}.\n"
            f"Expected: {reference_ids}\nObserved: {other_ids}"
        )


def parse_clipkit_keep_positions(path, expected_length):
    keep_positions = []
    for line_number, raw_line in enumerate(Path(path).read_text().splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 2:
            raise ValueError(f"Malformed ClipKIT log line {line_number}: {raw_line!r}")
        position = int(parts[0])
        if position < 1 or position > expected_length:
            raise ValueError(
                f"ClipKIT log position {position} is outside the protein alignment length {expected_length}."
            )
        if parts[1].lower() == "keep":
            keep_positions.append(position - 1)
    if not keep_positions:
        raise ValueError(f"ClipKIT log {path} did not report any retained columns.")
    if keep_positions != sorted(keep_positions):
        raise ValueError(f"ClipKIT log {path} is not ordered by ascending position.")
    return keep_positions


def derive_keep_positions_from_trimmed(original_records, trimmed_records):
    original_columns = list(zip(*(record["seq"] for record in original_records)))
    trimmed_columns = list(zip(*(record["seq"] for record in trimmed_records)))
    keep_positions = []
    trimmed_index = 0
    for original_index, column in enumerate(original_columns):
        if trimmed_index < len(trimmed_columns) and column == trimmed_columns[trimmed_index]:
            keep_positions.append(original_index)
            trimmed_index += 1
    if trimmed_index != len(trimmed_columns):
        raise ValueError(
            "Could not map all ClipKIT-retained protein columns back onto the original protein alignment."
        )
    return keep_positions


def slice_alignment_columns(records, keep_positions):
    return [
        {
            "id": record["id"],
            "description": record.get("description", record["id"]),
            "seq": "".join(record["seq"][index] for index in keep_positions),
        }
        for record in records
    ]


def slice_cds_by_protein_columns(records, keep_positions):
    trimmed_records = []
    for record in records:
        sequence = record["seq"]
        trimmed_records.append(
            {
                "id": record["id"],
                "description": record.get("description", record["id"]),
                "seq": "".join(sequence[index * 3:(index + 1) * 3] for index in keep_positions),
            }
        )
    return trimmed_records


def records_match(left_records, right_records):
    if len(left_records) != len(right_records):
        return False
    for left, right in zip(left_records, right_records):
        if left["id"] != right["id"] or left["seq"] != right["seq"]:
            return False
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--protein-alignment", required=True)
    parser.add_argument("--cds-alignment", required=True)
    parser.add_argument("--trimmed-protein-alignment", required=True)
    parser.add_argument("--protein-clipkit-log", default="")
    parser.add_argument("--direct-cds-clipkit-alignment", default="")
    parser.add_argument("--trimmed-cds-out", required=True)
    parser.add_argument("--report-out", required=True)
    parser.add_argument("--summary-out", required=True)
    args = parser.parse_args()

    protein_records = read_fasta_records(args.protein_alignment)
    cds_records = read_fasta_records(args.cds_alignment)
    trimmed_protein_records = read_fasta_records(args.trimmed_protein_alignment)

    protein_length = require_alignment(protein_records, "BABAPPAlign protein alignment")
    cds_length = require_alignment(cds_records, "BABAPPAlign CDS alignment")
    trimmed_protein_length = require_alignment(
        trimmed_protein_records, "ClipKIT-trimmed protein alignment"
    )

    require_same_ids(protein_records, cds_records, "BABAPPAlign protein alignment", "BABAPPAlign CDS alignment")
    require_same_ids(
        protein_records,
        trimmed_protein_records,
        "BABAPPAlign protein alignment",
        "ClipKIT-trimmed protein alignment",
    )

    if cds_length != protein_length * 3:
        raise ValueError(
            f"Aligned CDS length {cds_length} is not exactly three times the protein alignment length {protein_length}."
        )

    if args.protein_clipkit_log:
        keep_positions = parse_clipkit_keep_positions(args.protein_clipkit_log, protein_length)
        keep_source = "clipkit_log"
    else:
        keep_positions = derive_keep_positions_from_trimmed(protein_records, trimmed_protein_records)
        keep_source = "trimmed_alignment_subsequence"

    if len(keep_positions) != trimmed_protein_length:
        raise ValueError(
            f"Retained protein column count {len(keep_positions)} does not match the trimmed protein alignment length {trimmed_protein_length}."
        )

    reconstructed_trimmed_proteins = slice_alignment_columns(protein_records, keep_positions)
    if not records_match(reconstructed_trimmed_proteins, trimmed_protein_records):
        raise ValueError(
            "The projected protein columns do not reconstruct the ClipKIT-trimmed protein alignment. "
            "This indicates an inconsistent ClipKIT log or unexpected alignment reordering."
        )

    trimmed_cds_records = slice_cds_by_protein_columns(cds_records, keep_positions)
    write_fasta_records(trimmed_cds_records, args.trimmed_cds_out)

    direct_cds_present = bool(args.direct_cds_clipkit_alignment)
    direct_cds_length = 0
    direct_cds_exact_match = None
    if direct_cds_present:
        direct_cds_records = read_fasta_records(args.direct_cds_clipkit_alignment)
        direct_cds_length = require_alignment(direct_cds_records, "Direct ClipKIT CDS alignment")
        require_same_ids(
            cds_records,
            direct_cds_records,
            "BABAPPAlign CDS alignment",
            "Direct ClipKIT CDS alignment",
        )
        direct_cds_exact_match = records_match(trimmed_cds_records, direct_cds_records)
    trimmed_cds_length = len(trimmed_cds_records[0]["seq"])

    summary = {
        "status": "ok",
        "protein_alignment": args.protein_alignment,
        "cds_alignment": args.cds_alignment,
        "trimmed_protein_alignment": args.trimmed_protein_alignment,
        "protein_clipkit_log": args.protein_clipkit_log,
        "direct_cds_clipkit_alignment": args.direct_cds_clipkit_alignment,
        "projected_trimmed_cds_alignment": args.trimmed_cds_out,
        "keep_source": keep_source,
        "protein_columns_original": protein_length,
        "protein_columns_retained": len(keep_positions),
        "protein_columns_trimmed": protein_length - len(keep_positions),
        "cds_length_nt_original": cds_length,
        "cds_length_nt_retained": trimmed_cds_length,
        "retained_protein_positions_1based": [position + 1 for position in keep_positions],
        "direct_cds_clipkit_present": direct_cds_present,
        "direct_cds_clipkit_length_nt": direct_cds_length,
        "direct_cds_exact_match_to_projected_trim": direct_cds_exact_match,
    }

    Path(args.summary_out).write_text(json.dumps(summary, indent=2))

    report_lines = [
        "Protein-guided ClipKIT CDS audit summary",
        "========================================",
        f"Protein alignment: {args.protein_alignment}",
        f"CDS alignment: {args.cds_alignment}",
        f"Trimmed protein alignment: {args.trimmed_protein_alignment}",
        f"Projection strategy: {keep_source}",
        f"Original protein columns: {protein_length}",
        f"Retained protein columns: {len(keep_positions)}",
        f"Trimmed protein columns: {protein_length - len(keep_positions)}",
        f"Original CDS length (nt): {cds_length}",
        f"Projected CDS length (nt): {trimmed_cds_length}",
        f"Projected CDS QC output: {args.trimmed_cds_out}",
    ]
    if direct_cds_present:
        report_lines.extend(
            [
                f"Direct ClipKIT CDS alignment: {args.direct_cds_clipkit_alignment}",
                f"Direct ClipKIT CDS length (nt): {direct_cds_length}",
                "Direct ClipKIT CDS exact match to projected CDS trim: "
                f"{'yes' if direct_cds_exact_match else 'no'}",
            ]
        )
    Path(args.report_out).write_text("\n".join(report_lines) + "\n")


if __name__ == "__main__":
    main()
