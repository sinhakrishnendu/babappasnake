#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path

from common import as_bool, has_internal_stop_codon, read_fasta_records


def require_records(path, label):
    records = read_fasta_records(path)
    if not records:
        raise ValueError(f"{label} is empty: {path}")
    length = len(records[0]["seq"])
    ids = [record["id"] for record in records]
    if len(set(ids)) != len(ids):
        raise ValueError(f"{label} contains duplicate sequence identifiers.")
    for record in records:
        if len(record["seq"]) != length:
            raise ValueError(f"{label} contains inconsistent alignment lengths.")
    return records, length


def require_member_ids(records, expected_ids, label):
    observed_ids = [record["id"] for record in records]
    if observed_ids != expected_ids:
        raise ValueError(
            f"{label} sequence order does not match the selected orthogroup members.\n"
            f"Expected: {expected_ids}\nObserved: {observed_ids}"
        )


def validate_codon_alignment(records, expected_length, label):
    for record in records:
        if len(record["seq"]) != expected_length:
            raise ValueError(f"{label} has inconsistent alignment length for {record['id']}.")
        ungapped = record["seq"].replace("-", "")
        if len(ungapped) % 3 != 0:
            raise ValueError(f"{label} ungapped CDS length is not a multiple of three for {record['id']}.")
        if has_internal_stop_codon(record["seq"]):
            raise ValueError(f"{label} contains an internal stop codon for {record['id']}.")


def mean_gap_fraction(records, length):
    return round(
        sum(record["seq"].count("-") / length for record in records) / len(records),
        4,
    )


def path_or_empty(path):
    return path if path else ""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--babappalign-protein", required=True)
    parser.add_argument("--babappalign-cds", required=True)
    parser.add_argument("--trimmed-protein", default="")
    parser.add_argument("--trimmed-cds", default="")
    parser.add_argument("--projected-trimmed-cds", default="")
    parser.add_argument("--projection-summary", default="")
    parser.add_argument("--members-tsv", required=True)
    parser.add_argument("--clipkit-enabled", required=True)
    parser.add_argument("--iqtree-input-type", required=True, choices=["protein", "cds"])
    parser.add_argument("--effective-iqtree-alignment", required=True)
    parser.add_argument("--effective-codon-alignment", required=True)
    parser.add_argument("--json-out", required=True)
    parser.add_argument("--report-out", required=True)
    args = parser.parse_args()

    clipkit_enabled = as_bool(args.clipkit_enabled)

    with open(args.members_tsv) as handle:
        expected_ids = [row["member_id"] for row in csv.DictReader(handle, delimiter="\t")]

    bab_protein_records, bab_protein_length = require_records(
        args.babappalign_protein, "BABAPPAlign protein alignment"
    )
    bab_cds_records, bab_cds_length = require_records(args.babappalign_cds, "BABAPPAlign CDS alignment")
    require_member_ids(bab_protein_records, expected_ids, "BABAPPAlign protein alignment")
    require_member_ids(bab_cds_records, expected_ids, "BABAPPAlign CDS alignment")
    if bab_cds_length != bab_protein_length * 3:
        raise ValueError(
            f"BABAPPAlign CDS length {bab_cds_length} does not equal three times the BABAPPAlign protein length {bab_protein_length}."
        )
    validate_codon_alignment(bab_cds_records, bab_cds_length, "BABAPPAlign CDS alignment")

    trimmed_protein_records = None
    trimmed_protein_length = None
    trimmed_cds_records = None
    trimmed_cds_length = None
    projected_trimmed_cds_records = None
    projected_trimmed_cds_length = None
    projection_summary = {}

    if clipkit_enabled:
        if not args.trimmed_protein or not args.trimmed_cds:
            raise ValueError("ClipKIT is enabled but trimmed alignment paths were not provided.")

        trimmed_protein_records, trimmed_protein_length = require_records(
            args.trimmed_protein, "ClipKIT-trimmed protein alignment"
        )
        trimmed_cds_records, trimmed_cds_length = require_records(
            args.trimmed_cds, "Direct ClipKIT codon-trimmed CDS alignment"
        )
        require_member_ids(trimmed_protein_records, expected_ids, "ClipKIT-trimmed protein alignment")
        require_member_ids(trimmed_cds_records, expected_ids, "Direct ClipKIT codon-trimmed CDS alignment")

        if trimmed_protein_length > bab_protein_length:
            raise ValueError("Trimmed protein alignment is longer than the original BABAPPAlign protein alignment.")
        if trimmed_cds_length > bab_cds_length:
            raise ValueError("Trimmed CDS alignment is longer than the original BABAPPAlign CDS alignment.")
        validate_codon_alignment(
            trimmed_cds_records,
            trimmed_cds_length,
            "Direct ClipKIT codon-trimmed CDS alignment",
        )

        if args.projected_trimmed_cds:
            projected_trimmed_cds_records, projected_trimmed_cds_length = require_records(
                args.projected_trimmed_cds,
                "Protein-projected ClipKIT CDS alignment",
            )
            require_member_ids(
                projected_trimmed_cds_records,
                expected_ids,
                "Protein-projected ClipKIT CDS alignment",
            )
            validate_codon_alignment(
                projected_trimmed_cds_records,
                projected_trimmed_cds_length,
                "Protein-projected ClipKIT CDS alignment",
            )
            if projected_trimmed_cds_length > bab_cds_length:
                raise ValueError("Projected trimmed CDS alignment is longer than the original BABAPPAlign CDS alignment.")

        if args.projection_summary:
            projection_summary = json.loads(Path(args.projection_summary).read_text())
            expected_projected_length = projection_summary.get("cds_length_nt_retained")
            if expected_projected_length and expected_projected_length != projected_trimmed_cds_length:
                raise ValueError(
                    f"Projection summary reports CDS length {expected_projected_length}, but the projected CDS alignment is {projected_trimmed_cds_length} nt."
                )

    effective_codon_records, effective_codon_length = require_records(
        args.effective_codon_alignment,
        "Effective downstream CDS alignment",
    )
    require_member_ids(effective_codon_records, expected_ids, "Effective downstream CDS alignment")
    validate_codon_alignment(
        effective_codon_records,
        effective_codon_length,
        "Effective downstream CDS alignment",
    )

    effective_iqtree_records, effective_iqtree_length = require_records(
        args.effective_iqtree_alignment,
        "Effective IQ-TREE alignment",
    )
    require_member_ids(effective_iqtree_records, expected_ids, "Effective IQ-TREE alignment")

    summary = {
        "status": "ok",
        "clipkit_enabled": clipkit_enabled,
        "iqtree_input_type": args.iqtree_input_type,
        "effective_iqtree_alignment": args.effective_iqtree_alignment,
        "effective_codon_alignment": args.effective_codon_alignment,
        "downstream_cds_trimmed": clipkit_enabled,
        "sequence_count": len(expected_ids),
        "alignment_length_codons": effective_codon_length // 3,
        "babappalign": {
            "protein_alignment": args.babappalign_protein,
            "cds_alignment": args.babappalign_cds,
            "protein_length_aa": bab_protein_length,
            "cds_length_nt": bab_cds_length,
            "mean_protein_gap_fraction": mean_gap_fraction(bab_protein_records, bab_protein_length),
            "mean_cds_gap_fraction": mean_gap_fraction(bab_cds_records, bab_cds_length),
        },
        "clipkit": {
            "trimmed_protein_alignment": path_or_empty(args.trimmed_protein),
            "trimmed_cds_alignment": path_or_empty(args.trimmed_cds),
            "projected_trimmed_cds_alignment": path_or_empty(args.projected_trimmed_cds),
            "trimmed_protein_length_aa": trimmed_protein_length,
            "trimmed_cds_length_nt": trimmed_cds_length,
            "projected_trimmed_cds_length_nt": projected_trimmed_cds_length,
            "projection_summary": projection_summary,
            "cds_trim_strategy": "direct_clipkit_codon_mode",
        }
        if clipkit_enabled
        else {
            "trimmed_protein_alignment": "",
            "trimmed_cds_alignment": "",
            "projected_trimmed_cds_alignment": "",
            "trimmed_protein_length_aa": None,
            "trimmed_cds_length_nt": None,
            "projected_trimmed_cds_length_nt": None,
            "projection_summary": {},
            "cds_trim_strategy": "untrimmed",
        },
    }

    Path(args.json_out).write_text(json.dumps(summary, indent=2))

    report_lines = [
        "Alignment validation report",
        "=========================",
        f"ClipKIT enabled: {'yes' if clipkit_enabled else 'no'}",
        f"IQ-TREE input type: {args.iqtree_input_type}",
        f"BABAPPAlign protein alignment: {args.babappalign_protein}",
        f"BABAPPAlign CDS alignment: {args.babappalign_cds}",
        f"Effective IQ-TREE alignment: {args.effective_iqtree_alignment}",
        f"Effective downstream CDS alignment: {args.effective_codon_alignment}",
        f"Sequence count: {len(expected_ids)}",
        f"Effective downstream CDS length (codons): {effective_codon_length // 3}",
        "Validation: passed",
    ]
    if clipkit_enabled:
        report_lines.extend(
            [
                f"ClipKIT-trimmed protein alignment: {args.trimmed_protein}",
                f"ClipKIT codon-safe trimmed CDS alignment: {args.trimmed_cds}",
                "Primary CDS trim strategy: direct ClipKIT codon mode",
                f"Protein-projected CDS QC alignment: {args.projected_trimmed_cds or 'not generated'}",
            ]
        )
    Path(args.report_out).write_text("\n".join(report_lines) + "\n")


if __name__ == "__main__":
    main()
