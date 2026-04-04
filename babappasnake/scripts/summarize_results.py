#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def read_tsv(path: str) -> list[dict]:
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return []
    with open(p, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def pick_value(row: dict, keys: tuple[str, ...], default: str = "NA") -> str:
    for key in keys:
        value = row.get(key)
        if value is not None and str(value) != "":
            return str(value)
    return default


def detect_branchsite_significance_key(rows: list[dict]) -> str | None:
    preferred = ("significant_BH_0.05", "significant_BH_0.1", "significant_BH")
    available = {k for row in rows for k in row.keys()}
    for key in preferred:
        if key in available:
            return key
    for key in sorted(available):
        if key.startswith("significant_BH_"):
            return key
    return None


def count_meme_sites(meme_json_path: str, pcut: float) -> int:
    p = Path(meme_json_path)
    if not p.exists():
        return 0
    data = json.loads(p.read_text(encoding="utf-8"))
    content = json.dumps(data)
    return content.count('"p"')


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--rbh", required=True)
    p.add_argument("--mapping", required=True)
    p.add_argument("--absrel", required=True)
    p.add_argument("--absrel-meta", required=False)
    p.add_argument("--branchsite", required=True)
    p.add_argument("--hyphy-dir", required=True)
    p.add_argument("--gard-summary", required=False, default="")
    p.add_argument("--out", required=True)
    p.add_argument("--method", default="")
    p.add_argument("--trim-state", default="")
    p.add_argument("--meme-p", type=float, default=0.05)
    a = p.parse_args()

    rbh = read_tsv(a.rbh)
    mapping = read_tsv(a.mapping)
    absrel = read_tsv(a.absrel)
    branchsite = read_tsv(a.branchsite)
    hyphy_dir = Path(a.hyphy_dir)
    absrel_cutoff = 0.05
    if a.absrel_meta:
        mp = Path(a.absrel_meta)
        if mp.exists() and mp.stat().st_size > 0:
            meta = json.loads(mp.read_text(encoding="utf-8"))
            absrel_cutoff = float(meta.get("selected_pcut", absrel_cutoff))

    lines = []
    header = "BABAPPASNAKE EPISODIC SELECTION SUMMARY"
    if a.method:
        header += f" [{a.method}]"
    if a.trim_state:
        header += f" [{a.trim_state}]"
    lines.append(header)
    lines.append("=" * 42)
    lines.append("")
    lines.append(f"Orthogroup members recovered by RBH: {sum(1 for x in rbh if x.get('ortholog') not in {None, '', 'NA'}) + 1}")
    lines.append(f"CDS-to-protein mappings retained: {len(mapping)}")
    lines.append("")
    lines.append("Exploratory HyPhy stage")
    lines.append("-" * 24)
    if absrel:
        lines.append(f"Significant leaf branches by aBSREL (p <= {absrel_cutoff:.3f}): {len(absrel)}")
        for row in absrel:
            lines.append(f"  - {row['foreground_branch']} (aBSREL p={row['absrel_p']})")
    else:
        lines.append(f"No significant leaf branches detected by aBSREL at p <= {absrel_cutoff:.3f}.")
    lines.append("")
    lines.append("Branch-site codeml confirmation")
    lines.append("-" * 31)
    if branchsite:
        signif_key = detect_branchsite_significance_key(branchsite)
        for row in branchsite:
            significant_value = row.get(signif_key, "NA") if signif_key else "NA"
            lines.append(
                f"  - {pick_value(row, ('foreground_branch',), default='NA')}: "
                f"LRT={pick_value(row, ('LRT',), default='NA')}, "
                f"p={pick_value(row, ('p_value',), default='NA')}, "
                f"q={pick_value(row, ('q_value',), default='NA')}, "
                f"significant={significant_value}"
            )
    else:
        lines.append("No branch-site codeml tests were run because no foreground branch passed exploratory screening.")
    lines.append("")
    lines.append("Key files")
    lines.append("-" * 9)
    lines.append(f"RBH summary: {Path(a.rbh).resolve()}")
    lines.append(f"CDS mapping: {Path(a.mapping).resolve()}")
    lines.append(f"aBSREL foregrounds: {Path(a.absrel).resolve()}")
    lines.append(f"Branch-site results: {Path(a.branchsite).resolve()}")
    lines.append(f"HyPhy directory: {hyphy_dir.resolve()}")
    lines.append(f"Method: {a.method or 'NA'}")
    lines.append(f"Trim state: {a.trim_state or 'NA'}")
    lines.append(f"MEME threshold: {a.meme_p}")
    lines.append("")
    lines.append("Optional recombination screening (HyPhy GARD)")
    lines.append("-" * 44)
    gard_path = Path(str(a.gard_summary).strip()) if str(a.gard_summary).strip() else None
    if gard_path and gard_path.exists() and gard_path.stat().st_size > 0:
        try:
            gard = json.loads(gard_path.read_text(encoding="utf-8"))
            lines.append(f"GARD status: {gard.get('status', 'NA')}")
            lines.append(f"GARD breakpoints detected: {gard.get('breakpoints_detected', 'NA')}")
            lines.append(f"GARD breakpoint count: {gard.get('n_breakpoints', 'NA')}")
            note = str(gard.get("note", "")).strip()
            if note:
                lines.append(f"Note: {note}")
        except Exception:
            lines.append("GARD summary exists but could not be parsed.")
    else:
        lines.append("GARD not run for this pathway (default workflow behavior).")
    lines.append(
        "Caution: branch-site and downstream statistics reflect the full-length default pathway unless explicit fragment-aware routing is implemented."
    )
    Path(a.out).write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
