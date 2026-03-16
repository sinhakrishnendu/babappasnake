#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path

from common import as_bool, canonicalize_taxon, read_fasta_records, retention_cutoff, write_fasta_records


def choose_threshold(rows, min_ratio, max_loss):
    max_count = max(int(row["ortholog_count"]) for row in rows)
    cutoff = retention_cutoff(max_count, min_ratio, max_loss)
    acceptable = [row for row in rows if int(row["ortholog_count"]) >= cutoff]
    if acceptable:
        return max(acceptable, key=lambda row: int(row["threshold"])), cutoff, max_count
    best = max(rows, key=lambda row: (int(row["ortholog_count"]), int(row["threshold"])))
    return best, cutoff, max_count


def collapse_identical_members(members, fasta_records, outgroup_name):
    members_by_id = {row["member_id"]: row for row in members}
    order_index = {row["member_id"]: index for index, row in enumerate(members)}
    outgroup_key = canonicalize_taxon(outgroup_name)
    records_by_id = {record["id"]: record for record in fasta_records}

    sequence_groups = {}
    for record in fasta_records:
        sequence_groups.setdefault(record["seq"], []).append(record["id"])

    kept_ids = []
    duplicate_rows = []
    for member_ids in sequence_groups.values():
        if len(member_ids) == 1:
            kept_ids.append(member_ids[0])
            continue

        preferred_id = None
        for member_id in member_ids:
            row = members_by_id[member_id]
            if canonicalize_taxon(row["taxon_name"]) == outgroup_key:
                preferred_id = member_id
                break
        if preferred_id is None:
            preferred_id = min(member_ids, key=lambda member_id: order_index[member_id])
        kept_ids.append(preferred_id)
        for member_id in member_ids:
            if member_id == preferred_id:
                continue
            duplicate_rows.append(
                {
                    "kept_member_id": preferred_id,
                    "kept_taxon_name": members_by_id[preferred_id]["taxon_name"],
                    "removed_member_id": member_id,
                    "removed_taxon_name": members_by_id[member_id]["taxon_name"],
                }
            )

    kept_id_set = set(kept_ids)
    collapsed_members = [row for row in members if row["member_id"] in kept_id_set]
    collapsed_records = [records_by_id[row["member_id"]] for row in collapsed_members]
    return collapsed_members, collapsed_records, duplicate_rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary-tsv", required=True)
    parser.add_argument("--threshold-base-dir", required=True)
    parser.add_argument("--selected-dir", required=True)
    parser.add_argument("--query-meta", required=True)
    parser.add_argument("--outgroup-name", required=True)
    parser.add_argument("--min-ratio", required=True, type=float)
    parser.add_argument("--max-loss", required=True, type=int)
    parser.add_argument("--collapse-identical-proteins", action="store_true")
    parser.add_argument("--metadata-out", required=True)
    parser.add_argument("--report-out", required=True)
    parser.add_argument("--selected-members-out", required=True)
    parser.add_argument("--selected-fasta-out", required=True)
    args = parser.parse_args()

    with open(args.summary_tsv) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError("Threshold summary is empty; no orthogroup candidates were generated.")

    recommended, recommended_cutoff, recommended_max = choose_threshold(
        rows, args.min_ratio, args.max_loss
    )
    outgroup_rows = [row for row in rows if as_bool(row["outgroup_retained"])]
    if not outgroup_rows:
        raise ValueError(
            f"No tested coverage threshold retained the required outgroup {args.outgroup_name!r}. "
            "Review the threshold comparison outputs and proteome/outgroup naming."
        )

    selected, selected_cutoff, selected_max = choose_threshold(
        outgroup_rows, args.min_ratio, args.max_loss
    )
    selected_dir = Path(args.selected_dir)
    selected_dir.mkdir(parents=True, exist_ok=True)
    threshold_dir = Path(args.threshold_base_dir) / f"coverage_{selected['threshold']}"

    selection_reason = (
        f"Selected {selected['threshold']}% coverage because it is the highest threshold retaining "
        f"the outgroup while keeping at least {selected_cutoff}/{selected_max} orthologs "
        f"({args.min_ratio:.0%} of the outgroup-retaining maximum, maximum loss {args.max_loss})."
    )
    if selected["threshold"] != recommended["threshold"]:
        selection_reason += (
            f" The unconstrained recommendation was {recommended['threshold']}%, but it was not used "
            f"because the final orthogroup must retain the outgroup {args.outgroup_name}."
        )

    with open(threshold_dir / "orthogroup_members.tsv") as handle:
        members = list(csv.DictReader(handle, delimiter="\t"))
    fasta_records = read_fasta_records(threshold_dir / "orthogroup_proteins.faa")

    duplicate_rows = []
    if args.collapse_identical_proteins:
        members, fasta_records, duplicate_rows = collapse_identical_members(
            members, fasta_records, args.outgroup_name
        )

    with open(args.selected_members_out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(members[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(members)
    write_fasta_records(fasta_records, args.selected_fasta_out)

    metadata = {
        "orthogroup_id": members[0]["orthogroup_id"],
        "query_id": json.loads(Path(args.query_meta).read_text())["query_id"],
        "outgroup_name": args.outgroup_name,
        "selected_threshold": int(selected["threshold"]),
        "recommended_threshold_without_outgroup_constraint": int(recommended["threshold"]),
        "selection_reason": selection_reason,
        "selected_ortholog_count": int(selected["ortholog_count"]),
        "selected_member_count": len(members),
        "selected_member_count_before_collapse": int(selected["member_count"]),
        "collapse_identical_proteins": bool(args.collapse_identical_proteins),
        "collapsed_duplicate_count": len(duplicate_rows),
        "collapsed_duplicates": duplicate_rows,
        "retained_taxa": [row["taxon_name"] for row in members if row["member_type"] == "ortholog"],
        "outgroup_present": any(
            canonicalize_taxon(row["taxon_name"]) == canonicalize_taxon(args.outgroup_name)
            for row in members
        ),
    }
    Path(args.metadata_out).write_text(json.dumps(metadata, indent=2))

    report_lines = [
        f"Orthogroup ID: {metadata['orthogroup_id']}",
        f"Query ID: {metadata['query_id']}",
        f"Required outgroup: {args.outgroup_name}",
        f"Selected coverage threshold: {selected['threshold']}%",
        f"Automatic selection rationale: {selection_reason}",
        f"Orthologs retained (excluding query): {selected['ortholog_count']}",
        f"Members carried into downstream analysis: {len(members)}",
        (
            f"Collapsed identical protein sequences before CDS checkpoint: yes "
            f"({len(duplicate_rows)} redundant member(s) removed)"
            if duplicate_rows
            else "Collapsed identical protein sequences before CDS checkpoint: no"
        ),
        f"Outgroup retained: {'yes' if metadata['outgroup_present'] else 'no'}",
        "",
        "Orthogroup members:",
    ]
    for row in members:
        report_lines.append(
            f"  - {row.get('display_id', row['member_id'])} [{row['member_type']}] :: "
            f"{row['taxon_name']} :: {row['member_id']}"
        )
    if duplicate_rows:
        report_lines.extend(["", "Collapsed duplicate members:"])
        for row in duplicate_rows:
            report_lines.append(
                f"  - kept {row['kept_member_id']} ({row['kept_taxon_name']}); "
                f"removed {row['removed_member_id']} ({row['removed_taxon_name']})"
            )
    Path(args.report_out).write_text("\n".join(report_lines) + "\n")


if __name__ == "__main__":
    main()
