#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path


def parse_absrel(data):
    tested = data.get("tested", {}).get("0", {})
    branch_attributes = data.get("branch attributes", {}).get("0", {})
    threshold = data.get("test results", {}).get("P-value threshold", 0.05)
    rows = []
    for branch_name, attrs in branch_attributes.items():
        if branch_name == "attributes":
            continue
        role = tested.get(branch_name, "untested")
        p_uncorrected = attrs.get("Uncorrected P-value", attrs.get("p-value"))
        p_corrected = attrs.get("Corrected P-value")
        significant_raw = role == "test" and p_uncorrected is not None and p_uncorrected <= threshold
        significant_corrected = (
            role == "test"
            and (
                (p_corrected is not None and p_corrected <= threshold)
                or (p_corrected is None and significant_raw)
            )
        )
        rows.append(
            {
                "branch_name": branch_name,
                "original_name": attrs.get("original name", branch_name),
                "tested_role": role,
                "lrt": attrs.get("LRT", ""),
                "p_uncorrected": p_uncorrected if p_uncorrected is not None else "",
                "p_corrected": p_corrected if p_corrected is not None else "",
                "significant_raw": "yes" if significant_raw else "no",
                "significant_corrected": "yes" if significant_corrected else "no",
            }
        )
    rows.sort(key=lambda row: (row["tested_role"] != "test", float(row["p_uncorrected"] or 1.0)))
    summary = {
        "positive_test_results": data.get("test results", {}).get("positive test results", 0),
        "tested_branches": data.get("test results", {}).get("tested", 0),
        "pvalue_threshold": threshold,
    }
    return rows, summary


def parse_busted(data):
    results = data.get("test results", {})
    pvalue = results.get("p-value")
    summary = {
        "lrt": results.get("LRT"),
        "p_value": pvalue,
        "significant": pvalue is not None and pvalue <= 0.05,
    }
    return [], summary


def normalize_meme_header(header):
    header = header.replace("&alpha;", "alpha")
    header = header.replace("&beta;<sup>1</sup>", "beta1")
    header = header.replace("p<sup>1</sup>", "p1")
    header = header.replace("&beta;<sup>+</sup>", "beta_plus")
    header = header.replace("p<sup>+</sup>", "p_plus")
    header = header.replace("# branches under selection", "branches_under_selection")
    header = header.replace("Total branch length", "total_branch_length")
    header = header.replace("MEME LogL", "meme_logl")
    header = header.replace("FEL LogL", "fel_logl")
    header = header.replace("FEL &alpha;", "fel_alpha")
    header = header.replace("FEL &beta;", "fel_beta")
    header = header.replace("p-value", "p_value")
    header = header.replace(" ", "_").replace("-", "_").lower()
    return header


def parse_meme(data, significance_threshold):
    headers = data["MLE"]["headers"]
    header_names = [normalize_meme_header(item[0] if isinstance(item, list) else item) for item in headers]
    content = data["MLE"]["content"]
    rows = []
    if isinstance(content, dict):
        partitions = content.values()
    else:
        partitions = [content]
    site_index = 1
    for partition_rows in partitions:
        for values in partition_rows:
            row = dict(zip(header_names, values))
            row["site"] = site_index
            row["significant"] = "yes" if row.get("p_value", 1.0) <= significance_threshold else "no"
            rows.append(row)
            site_index += 1
    summary = {
        "significant_sites": sum(1 for row in rows if row["significant"] == "yes"),
        "site_count": len(rows),
        "pvalue_threshold": significance_threshold,
    }
    return rows, summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", choices=["absrel", "busted", "meme"], required=True)
    parser.add_argument("--json-in", required=True)
    parser.add_argument("--table-out", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--meme-threshold", type=float, default=0.1)
    args = parser.parse_args()

    data = json.loads(Path(args.json_in).read_text())
    if args.method == "absrel":
        rows, summary = parse_absrel(data)
    elif args.method == "busted":
        rows, summary = parse_busted(data)
    else:
        rows, summary = parse_meme(data, args.meme_threshold)

    Path(args.summary_out).write_text(json.dumps(summary, indent=2))

    if rows:
        with open(args.table_out, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    else:
        Path(args.table_out).write_text("")


if __name__ == "__main__":
    main()
