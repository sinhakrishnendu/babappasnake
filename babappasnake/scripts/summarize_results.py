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
    p.add_argument("--out", required=True)
    p.add_argument("--meme-p", type=float, default=0.1)
    a = p.parse_args()

    rbh = read_tsv(a.rbh)
    mapping = read_tsv(a.mapping)
    absrel = read_tsv(a.absrel)
    branchsite = read_tsv(a.branchsite)
    hyphy_dir = Path(a.hyphy_dir)
    absrel_cutoff = 0.1
    if a.absrel_meta:
        mp = Path(a.absrel_meta)
        if mp.exists() and mp.stat().st_size > 0:
            meta = json.loads(mp.read_text(encoding="utf-8"))
            absrel_cutoff = float(meta.get("selected_pcut", absrel_cutoff))

    lines = []
    lines.append("BABAPPASNAKE EPISODIC SELECTION SUMMARY")
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
        for row in branchsite:
            lines.append(
                f"  - {row['foreground_branch']}: LRT={row['LRT']}, p={row['p_value']}, q={row['q_value']}, significant={row['significant_BH_0.1']}"
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
    Path(a.out).write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
