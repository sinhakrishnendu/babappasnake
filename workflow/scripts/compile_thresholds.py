#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path

from common import (
    canonicalize_taxon,
    make_unique_member_ids,
    read_fasta_records,
    write_fasta_records,
)


def load_manifest(path):
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def load_candidate_rows(paths):
    rows_by_species = {}
    for path in paths:
        with open(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            rows = []
            for row in reader:
                row["rank_by_bitscore"] = int(row["rank_by_bitscore"])
                row["bitscore"] = float(row["bitscore"])
                row["evalue"] = float(row["evalue"])
                row["pident"] = float(row["pident"])
                row["alignment_length"] = int(float(row["alignment_length"]))
                row["query_length"] = int(float(row["query_length"]))
                row["subject_length"] = int(float(row["subject_length"]))
                row["query_coverage"] = float(row["query_coverage"])
                row["subject_coverage"] = float(row["subject_coverage"])
                row["reverse_bitscore"] = float(row["reverse_bitscore"])
                rows.append(row)
        rows.sort(key=lambda item: (item["rank_by_bitscore"], -item["bitscore"]))
        species_id = Path(path).parent.name
        rows_by_species[species_id] = rows
    return rows_by_species


def extract_protein_sequence(proteome_path, record_id):
    for record in read_fasta_records(proteome_path):
        if record["id"] == record_id:
            return record["seq"]
    raise KeyError(f"Sequence {record_id} not found in {proteome_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--query-fasta", required=True)
    parser.add_argument("--query-meta", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--outgroup-name", required=True)
    parser.add_argument("--thresholds-json", required=True)
    parser.add_argument("--candidate-files", nargs="+", required=True)
    parser.add_argument("--threshold-base-dir", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--summary-json-out", required=True)
    args = parser.parse_args()

    thresholds = json.loads(args.thresholds_json)
    query_meta = json.loads(Path(args.query_meta).read_text())
    query_record = read_fasta_records(args.query_fasta)[0]
    manifest_rows = load_manifest(args.manifest)
    manifest_by_species = {row["species_id"]: row for row in manifest_rows}
    candidate_rows = load_candidate_rows(args.candidate_files)
    manifest_taxa = [row["taxon_name"] for row in manifest_rows]
    outgroup_key = canonicalize_taxon(args.outgroup_name)
    orthogroup_id = f"{query_meta['query_member_id']}_rbh_orthogroup"

    summary_rows = []
    threshold_base_dir = Path(args.threshold_base_dir)
    threshold_base_dir.mkdir(parents=True, exist_ok=True)

    for threshold in thresholds:
        threshold_fraction = float(threshold) / 100.0
        selected_rows = []
        for species_id, manifest_row in manifest_by_species.items():
            picked = None
            for row in candidate_rows.get(species_id, []):
                if (
                    row["query_coverage"] >= threshold_fraction
                    and row["subject_coverage"] >= threshold_fraction
                ):
                    picked = row.copy()
                    break
            if picked:
                selected_rows.append(picked)

        member_labels = [row["taxon_name"] for row in selected_rows]
        unique_member_ids = make_unique_member_ids(member_labels)

        selected_dir = threshold_base_dir / f"coverage_{int(threshold)}"
        selected_dir.mkdir(parents=True, exist_ok=True)

        member_rows = []
        fasta_records = []

        retained_taxa = []
        for member_id, row in zip(unique_member_ids, selected_rows):
            retained_taxa.append(row["taxon_name"])
            original_id = row["subject_id"]
            member_rows.append(
                {
                    "orthogroup_id": orthogroup_id,
                    "member_id": original_id,
                    "display_id": member_id,
                    "member_type": "ortholog",
                    "species_id": row["species_id"],
                    "taxon_name": row["taxon_name"],
                    "original_id": original_id,
                    "query_coverage": f"{row['query_coverage']:.4f}",
                    "subject_coverage": f"{row['subject_coverage']:.4f}",
                    "bitscore": f"{row['bitscore']:.4f}",
                    "pident": f"{row['pident']:.4f}",
                    "threshold": threshold,
                    "is_outgroup": "yes"
                    if canonicalize_taxon(row["taxon_name"]) == outgroup_key
                    else "no",
                }
            )
            fasta_records.append(
                {
                    "id": original_id,
                    "seq": extract_protein_sequence(row["proteome_fasta"], original_id),
                }
            )

        with open(selected_dir / "orthogroup_members.tsv", "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(member_rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(member_rows)
        write_fasta_records(fasta_records, selected_dir / "orthogroup_proteins.faa")

        missing_taxa = [
            taxon_name for taxon_name in manifest_taxa if taxon_name not in retained_taxa
        ]
        outgroup_retained = any(
            canonicalize_taxon(taxon_name) == outgroup_key for taxon_name in retained_taxa
        )
        summary_rows.append(
            {
                "threshold": threshold,
                "orthogroup_id": orthogroup_id,
                "ortholog_count": len(selected_rows),
                "member_count": len(member_rows),
                "outgroup_retained": str(outgroup_retained).lower(),
                "retained_taxa": "; ".join(retained_taxa),
                "missing_taxa": "; ".join(missing_taxa),
            }
        )

    with open(args.summary_out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    Path(args.summary_json_out).write_text(
        json.dumps(
            {
                "query_id": query_meta["query_id"],
                "query_length_aa": len(query_record["seq"]),
                "orthogroup_id": orthogroup_id,
                "thresholds": summary_rows,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
