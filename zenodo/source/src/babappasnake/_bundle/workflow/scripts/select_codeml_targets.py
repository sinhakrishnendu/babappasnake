#!/usr/bin/env python3
import argparse
import csv
import json
import re
from pathlib import Path


def is_leaf_branch(row):
    return not re.fullmatch(r"Node\d+", row["branch_name"] or "")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--absrel-table", required=True)
    parser.add_argument("--fallback-count", required=True, type=int)
    parser.add_argument("--selection-out", required=True)
    parser.add_argument("--summary-out", required=True)
    args = parser.parse_args()

    with open(args.absrel_table) as handle:
        rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("tested_role") == "test" and is_leaf_branch(row)
        ]
    if not rows:
        raise ValueError("No testable leaf branches were present in the aBSREL summary.")

    for row in rows:
        row["p_uncorrected_float"] = float(row["p_uncorrected"] or 1.0)
        row["p_corrected_float"] = float(row["p_corrected"] or row["p_uncorrected"] or 1.0)

    significant = [row for row in rows if row["significant_corrected"] == "yes"]
    selected = []
    if significant:
        significant.sort(key=lambda row: row["p_corrected_float"])
        for rank, row in enumerate(significant, start=1):
            selected.append(
                {
                    "branch_id": row["branch_name"],
                    "branch_name": row["original_name"] or row["branch_name"],
                    "original_name": row["original_name"],
                    "selection_strategy": "significant_absrel",
                    "absrel_rank": rank,
                    "absrel_p_uncorrected": row["p_uncorrected"],
                    "absrel_p_corrected": row["p_corrected"],
                    "absrel_lrt": row["lrt"],
                    "reason": "aBSREL significant after branch-wise correction",
                }
            )
        strategy_summary = (
            "Used all aBSREL-significant terminal branches for codeml foreground testing; "
            "internal nodes were excluded."
        )
    else:
        fallback = sorted(rows, key=lambda row: row["p_uncorrected_float"])[: args.fallback_count]
        for rank, row in enumerate(fallback, start=1):
            selected.append(
                {
                    "branch_id": row["branch_name"],
                    "branch_name": row["original_name"] or row["branch_name"],
                    "original_name": row["original_name"],
                    "selection_strategy": "fallback_lowest_absrel_p",
                    "absrel_rank": rank,
                    "absrel_p_uncorrected": row["p_uncorrected"],
                    "absrel_p_corrected": row["p_corrected"],
                    "absrel_lrt": row["lrt"],
                    "reason": "No significant aBSREL branches; selected among the lowest uncorrected p-values.",
                }
            )
        strategy_summary = (
            f"No significant aBSREL branches were detected; selected the top {len(selected)} "
            "lowest-p terminal branches for branch-site follow-up, excluding internal nodes."
        )

    with open(args.selection_out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(selected[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(selected)

    Path(args.summary_out).write_text(
        json.dumps(
            {
                "selection_strategy_summary": strategy_summary,
                "selected_branch_count": len(selected),
                "branches": selected,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
