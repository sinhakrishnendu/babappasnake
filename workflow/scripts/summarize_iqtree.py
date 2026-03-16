#!/usr/bin/env python3
import argparse
import json
import re
from pathlib import Path


def extract(pattern, text, group=1, default=""):
    match = re.search(pattern, text, re.MULTILINE)
    return match.group(group) if match else default


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--iqtree-report", required=True)
    parser.add_argument("--treefile", required=True)
    parser.add_argument("--alignment", required=True)
    parser.add_argument("--input-type", required=True, choices=["protein", "cds"])
    parser.add_argument("--sequence-type", required=True)
    parser.add_argument("--summary-out", required=True)
    args = parser.parse_args()

    text = Path(args.iqtree_report).read_text(errors="ignore")
    summary = {
        "treefile": args.treefile,
        "iqtree_report": args.iqtree_report,
        "alignment": args.alignment,
        "input_type": args.input_type,
        "sequence_type": args.sequence_type,
        "best_fit_model": extract(r"Best-fit model according to .*?:\s+([^\n]+)", text),
        "log_likelihood": extract(r"Log-likelihood of the tree:\s+([-\d.]+)", text),
        "note": extract(r"(NOTE:\s+[^\n]+)", text),
    }
    Path(args.summary_out).write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
