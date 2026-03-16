#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

from common import canonicalize_taxon, read_fasta_records, sanitize_identifier


def load_species_metadata(path):
    if not path:
        return {}
    import csv

    metadata = {}
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "taxon_name" not in reader.fieldnames:
            raise ValueError("species_metadata_tsv must contain a taxon_name column.")
        for row in reader:
            key = row.get("species_id") or sanitize_identifier(Path(row["proteome_path"]).stem)
            metadata[key] = row["taxon_name"]
    return metadata


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--query-fasta", required=True)
    parser.add_argument("--outgroup-name", required=True)
    parser.add_argument("--query-taxon-name", default="")
    parser.add_argument("--species-metadata-tsv", default="")
    parser.add_argument("--proteome-json", required=True)
    parser.add_argument("--manifest-out", required=True)
    parser.add_argument("--query-meta-out", required=True)
    parser.add_argument("--validation-out", required=True)
    args = parser.parse_args()

    proteomes = json.loads(args.proteome_json)
    species_metadata = load_species_metadata(args.species_metadata_tsv)

    query_records = read_fasta_records(args.query_fasta)
    if len(query_records) != 1:
        raise ValueError(
            f"Query FASTA must contain exactly one sequence, found {len(query_records)} in {args.query_fasta}"
        )
    query_record = query_records[0]
    query_member_id = f"query__{sanitize_identifier(query_record['id'])}"
    query_taxon_name = args.query_taxon_name or query_record["id"]

    manifest_path = Path(args.manifest_out)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    import csv

    candidate_taxa = []
    with open(manifest_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "species_id",
                "file_label",
                "proteome_path",
                "taxon_name",
                "taxon_key",
            ],
        )
        writer.writeheader()
        for record in proteomes:
            proteome_path = Path(record["path"])
            if not proteome_path.exists():
                raise FileNotFoundError(f"Proteome FASTA not found: {proteome_path}")
            first_record = read_fasta_records(proteome_path)
            if not first_record:
                raise ValueError(f"Proteome FASTA is empty: {proteome_path}")
            taxon_name = species_metadata.get(record["species_id"], record["file_label"])
            candidate_taxa.append(taxon_name)
            writer.writerow(
                {
                    "species_id": record["species_id"],
                    "file_label": record["file_label"],
                    "proteome_path": str(proteome_path),
                    "taxon_name": taxon_name,
                    "taxon_key": canonicalize_taxon(taxon_name),
                }
            )

    query_meta = {
        "query_id": query_record["id"],
        "query_member_id": query_member_id,
        "query_taxon_name": query_taxon_name,
        "query_length_aa": len(query_record["seq"]),
        "outgroup_name": args.outgroup_name,
    }
    Path(args.query_meta_out).write_text(json.dumps(query_meta, indent=2))

    outgroup_seen_in_manifest = canonicalize_taxon(args.outgroup_name) in {
        canonicalize_taxon(value) for value in candidate_taxa
    }
    validation = {
        "query_fasta": args.query_fasta,
        "query_id": query_record["id"],
        "proteome_count": len(proteomes),
        "outgroup_name": args.outgroup_name,
        "outgroup_seen_in_manifest": outgroup_seen_in_manifest,
        "note": (
            "Outgroup matching in downstream selection also uses RBH-derived taxon labels. "
            "If your filenames are not taxon names, use species_metadata_tsv for explicit mapping."
        ),
    }
    Path(args.validation_out).write_text(json.dumps(validation, indent=2))


if __name__ == "__main__":
    main()
