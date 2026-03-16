#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

from Bio import Phylo

from common import canonicalize_taxon


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--treefile", required=True)
    parser.add_argument("--members-tsv", required=True)
    parser.add_argument("--outgroup-name", required=True)
    parser.add_argument("--rooted-tree-out", required=True)
    parser.add_argument("--branch-map-out", required=True)
    args = parser.parse_args()

    with open(args.members_tsv) as handle:
        members = list(csv.DictReader(handle, delimiter="\t"))

    outgroup_key = canonicalize_taxon(args.outgroup_name)
    outgroup_member = None
    member_taxa = {}
    for row in members:
        member_taxa[row["member_id"]] = row["taxon_name"]
        if row["member_type"] == "ortholog" and canonicalize_taxon(row["taxon_name"]) == outgroup_key:
            outgroup_member = row["member_id"]
    if outgroup_member is None:
        raise ValueError(f"Outgroup {args.outgroup_name!r} is absent from the selected orthogroup.")

    tree = Phylo.read(args.treefile, "newick")
    matches = [clade for clade in tree.find_clades() if clade.name == outgroup_member]
    if not matches:
        raise ValueError(f"Outgroup member {outgroup_member!r} was not found in the IQ-TREE tree.")

    tree.root_with_outgroup(matches[0])
    counter = 1
    branch_rows = []
    for clade in tree.find_clades(order="preorder"):
        clade.confidence = None
        if clade.is_terminal():
            branch_rows.append(
                {
                    "branch_name": clade.name,
                    "branch_type": "leaf",
                    "taxon_name": member_taxa.get(clade.name, ""),
                }
            )
            continue
        original_name = clade.name or ""
        clade.name = f"Node{counter}"
        counter += 1
        branch_rows.append(
            {
                "branch_name": clade.name,
                "branch_type": "root" if clade == tree.root else "internal",
                "taxon_name": "",
                "original_label": original_name,
            }
        )

    Path(args.rooted_tree_out).parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, args.rooted_tree_out, "newick")

    with open(args.branch_map_out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(branch_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(branch_rows)


if __name__ == "__main__":
    main()
