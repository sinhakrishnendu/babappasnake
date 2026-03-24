#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from itertools import combinations
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def parse_methods(raw: str) -> list[str]:
    methods = [token.strip() for token in raw.split(",") if token.strip()]
    if not methods:
        raise RuntimeError("No alignment methods were provided.")
    deduped: list[str] = []
    seen = set()
    for method in methods:
        if method not in seen:
            seen.add(method)
            deduped.append(method)
    return deduped


def parse_bool(value: str) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def detect_significance_key(rows: list[dict[str, str]]) -> str | None:
    preferred = ("significant_BH_0.05", "significant_BH_0.1", "significant_BH")
    available = {key for row in rows for key in row.keys()}
    for key in preferred:
        if key in available:
            return key
    for key in sorted(available):
        if key.startswith("significant_BH_"):
            return key
    return None


def jaccard(a: set[str], b: set[str]) -> float:
    union = a | b
    if not union:
        return 1.0
    return len(a & b) / len(union)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", required=True)
    p.add_argument("--methods", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()

    outdir = Path(a.outdir)
    methods = parse_methods(a.methods)

    absrel_sets: dict[str, set[str]] = {}
    branchsite_sig_sets: dict[str, set[str]] = {}
    branchsite_all_sets: dict[str, set[str]] = {}
    absrel_cutoffs: dict[str, float | None] = {}

    for method in methods:
        absrel_tsv = outdir / "hyphy" / method / "significant_foregrounds.tsv"
        absrel_meta = outdir / "hyphy" / method / "foreground_threshold.json"
        branchsite_tsv = outdir / "branchsite" / method / "branchsite_results.tsv"

        absrel_rows = read_tsv(absrel_tsv)
        branchsite_rows = read_tsv(branchsite_tsv)
        signif_key = detect_significance_key(branchsite_rows)

        absrel_sets[method] = {
            str(row.get("foreground_branch", "")).strip()
            for row in absrel_rows
            if str(row.get("foreground_branch", "")).strip()
        }
        branchsite_all_sets[method] = {
            str(row.get("foreground_branch", "")).strip()
            for row in branchsite_rows
            if str(row.get("foreground_branch", "")).strip()
        }
        if signif_key:
            branchsite_sig_sets[method] = {
                str(row.get("foreground_branch", "")).strip()
                for row in branchsite_rows
                if str(row.get("foreground_branch", "")).strip() and parse_bool(row.get(signif_key, ""))
            }
        else:
            branchsite_sig_sets[method] = set()

        if absrel_meta.exists() and absrel_meta.stat().st_size > 0:
            meta = json.loads(absrel_meta.read_text(encoding="utf-8"))
            try:
                absrel_cutoffs[method] = float(meta.get("selected_pcut"))
            except Exception:
                absrel_cutoffs[method] = None
        else:
            absrel_cutoffs[method] = None

    lines: list[str] = []
    lines.append("BABAPPASNAKE ALIGNMENT METHOD REPRODUCIBILITY SUMMARY")
    lines.append("=" * 54)
    lines.append("")
    lines.append(f"Methods compared: {', '.join(methods)}")
    lines.append("")
    lines.append("Per-method hit counts")
    lines.append("-" * 21)
    for method in methods:
        cutoff = absrel_cutoffs.get(method)
        cutoff_text = f"{cutoff:.3f}" if cutoff is not None else "NA"
        lines.append(
            f"- {method}: aBSREL={len(absrel_sets[method])} (cutoff={cutoff_text}), "
            f"branch-site tested={len(branchsite_all_sets[method])}, "
            f"branch-site significant={len(branchsite_sig_sets[method])}"
        )

    lines.append("")
    lines.append("Pairwise overlap (Jaccard)")
    lines.append("-" * 26)
    for left, right in combinations(methods, 2):
        absrel_j = jaccard(absrel_sets[left], absrel_sets[right])
        branch_j = jaccard(branchsite_sig_sets[left], branchsite_sig_sets[right])
        lines.append(
            f"- {left} vs {right}: aBSREL={absrel_j:.3f}, branch-site-significant={branch_j:.3f}"
        )
    if len(methods) < 2:
        lines.append("- Only one method selected; pairwise comparison not applicable.")

    lines.append("")
    lines.append("Consensus and unique branches")
    lines.append("-" * 29)
    absrel_consensus = set.intersection(*(absrel_sets[m] for m in methods)) if methods else set()
    branch_consensus = set.intersection(*(branchsite_sig_sets[m] for m in methods)) if methods else set()
    lines.append(f"- Consensus aBSREL branches across all methods: {len(absrel_consensus)}")
    if absrel_consensus:
        for branch in sorted(absrel_consensus):
            lines.append(f"  - {branch}")
    lines.append(f"- Consensus branch-site significant branches across all methods: {len(branch_consensus)}")
    if branch_consensus:
        for branch in sorted(branch_consensus):
            lines.append(f"  - {branch}")

    for method in methods:
        other_absrel = set().union(*(absrel_sets[m] for m in methods if m != method))
        other_branch = set().union(*(branchsite_sig_sets[m] for m in methods if m != method))
        unique_absrel = absrel_sets[method] - other_absrel
        unique_branch = branchsite_sig_sets[method] - other_branch
        lines.append(
            f"- {method} unique: aBSREL={len(unique_absrel)}, branch-site-significant={len(unique_branch)}"
        )

    Path(a.out).write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
