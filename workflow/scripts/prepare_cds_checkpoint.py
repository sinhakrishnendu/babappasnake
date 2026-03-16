#!/usr/bin/env python3
import argparse
import csv
import json
import shutil
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--members-tsv", required=True)
    parser.add_argument("--proteins-fasta", required=True)
    parser.add_argument("--selection-metadata", required=True)
    parser.add_argument("--checkpoint-dir", required=True)
    parser.add_argument("--seed-cds-fasta", default="")
    args = parser.parse_args()

    checkpoint_dir = Path(args.checkpoint_dir)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(args.members_tsv, checkpoint_dir / "orthogroup_members.tsv")
    shutil.copyfile(args.proteins_fasta, checkpoint_dir / "orthogroup_proteins.faa")
    selection = json.loads(Path(args.selection_metadata).read_text())

    with open(args.members_tsv) as handle:
        members = list(csv.DictReader(handle, delimiter="\t"))
    required_headers = [row["member_id"] for row in members]
    expected_cds_path = checkpoint_dir / "user_supplied_cds.fasta"

    metadata = {
        "orthogroup_id": selection["orthogroup_id"],
        "selected_threshold": selection["selected_threshold"],
        "expected_cds_fasta": str(expected_cds_path),
        "required_member_ids": required_headers,
        "selection_reason": selection["selection_reason"],
    }
    (checkpoint_dir / "checkpoint_metadata.json").write_text(json.dumps(metadata, indent=2))
    (checkpoint_dir / "required_member_ids.txt").write_text("\n".join(required_headers) + "\n")
    (checkpoint_dir / "README.txt").write_text(
        "\n".join(
            [
                f"Orthogroup checkpoint: {selection['orthogroup_id']}",
                f"Selected threshold: {selection['selected_threshold']}%",
                "",
                "Place a CDS FASTA covering every listed member into:",
                f"  {expected_cds_path}",
                "",
                "Rules for the user-supplied CDS file:",
                "  1. Headers must exactly match the member_id values in orthogroup_members.tsv.",
                "  2. Supply one CDS sequence per downstream analysis member.",
                "  3. CDS records must be unaligned, in-frame, and free of internal stop codons.",
                "",
                "Resume by rerunning the same snakemake command after the file is present.",
            ]
        )
        + "\n"
    )

    if args.seed_cds_fasta:
        seed_path = Path(args.seed_cds_fasta)
        if seed_path.exists():
            shutil.copyfile(seed_path, expected_cds_path)


if __name__ == "__main__":
    main()
