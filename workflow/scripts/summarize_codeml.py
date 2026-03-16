#!/usr/bin/env python3
import argparse
import csv
import json
import re
from pathlib import Path

from common import benjamini_hochberg, mixed_chi_square_pvalue


def parse_lnl(path):
    text = Path(path).read_text(errors="ignore")
    match = re.search(r"lnL\(ntime:\s*\d+\s+np:\s*(\d+)\):\s*([-\d.]+)", text)
    if not match:
        raise ValueError(f"Could not parse lnL from {path}")
    return int(match.group(1)), float(match.group(2)), text


def parse_beb_sites(text):
    sites = []
    capture = False
    for line in text.splitlines():
        if "Bayes Empirical Bayes (BEB) analysis" in line:
            capture = False
            continue
        if "Positive sites for foreground lineages Prob(w>1):" in line:
            capture = True
            continue
        if capture:
            if not line.strip() or line.startswith("The grid"):
                break
            match = re.match(r"\s*(\d+)\s+([A-Z\*\-])\s+([0-9.]+)\*{0,2}", line)
            if match:
                sites.append(
                    {
                        "site": int(match.group(1)),
                        "residue": match.group(2),
                        "posterior_probability": float(match.group(3)),
                    }
                )
    return sites


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--selection-tsv", required=True)
    parser.add_argument("--codeml-dir", required=True)
    parser.add_argument("--fdr-alpha", required=True, type=float)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--summary-json-out", required=True)
    args = parser.parse_args()

    codeml_dir = Path(args.codeml_dir)
    with open(args.selection_tsv) as handle:
        selection_rows = list(csv.DictReader(handle, delimiter="\t"))
    summary_rows = []
    raw_pvalues = []
    for row in selection_rows:
        branch_dir = codeml_dir / row["branch_id"]
        np_alt, lnl_alt, alt_text = parse_lnl(branch_dir / "branch_site_alt" / "output.txt")
        np_null, lnl_null, null_text = parse_lnl(branch_dir / "branch_site_null" / "output.txt")
        lrt = max(0.0, 2.0 * (lnl_alt - lnl_null))
        raw_p = mixed_chi_square_pvalue(lrt)
        raw_pvalues.append(raw_p)
        beb_sites = parse_beb_sites(alt_text)
        summary_rows.append(
            {
                "branch_id": row["branch_id"],
                "branch_name": row["branch_name"],
                "selection_strategy": row["selection_strategy"],
                "absrel_rank": row["absrel_rank"],
                "absrel_p_uncorrected": row["absrel_p_uncorrected"],
                "absrel_p_corrected": row["absrel_p_corrected"],
                "lnL_alt": f"{lnl_alt:.6f}",
                "np_alt": np_alt,
                "lnL_null": f"{lnl_null:.6f}",
                "np_null": np_null,
                "lrt_statistic": f"{lrt:.6f}",
                "raw_p_value": raw_p,
                "beb_sites": "; ".join(
                    f"{site['site']}:{site['residue']}:{site['posterior_probability']:.3f}"
                    for site in beb_sites
                ),
                "beb_sites_ge_0_95": "; ".join(
                    str(site["site"]) for site in beb_sites if site["posterior_probability"] >= 0.95
                ),
            }
        )

    adjusted = benjamini_hochberg(raw_pvalues)
    for row, adjusted_p in zip(summary_rows, adjusted):
        row["bh_fdr"] = adjusted_p
        row["significant_after_bh"] = "yes" if adjusted_p <= args.fdr_alpha else "no"

    with open(args.summary_out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    Path(args.summary_json_out).write_text(
        json.dumps(
            {
                "fdr_alpha": args.fdr_alpha,
                "tested_branches": len(summary_rows),
                "significant_branches_after_bh": [
                    row["branch_name"] for row in summary_rows if row["significant_after_bh"] == "yes"
                ],
                "rows": summary_rows,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
