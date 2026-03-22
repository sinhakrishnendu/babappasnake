#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import re
from collections import defaultdict


def canonicalize_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")


def extract_tip_names(tree: str) -> list[str]:
    # Leaf labels are names that appear after "(" or "," and before ":".
    return re.findall(r"(?:(?<=\()|(?<=,))([^():,]+)(?=:)", tree)


def resolve_foreground_name(tree: str, foreground: str) -> str:
    tip_names = extract_tip_names(tree)
    if foreground in tip_names:
        return foreground

    by_canon: dict[str, list[str]] = defaultdict(list)
    for tip in tip_names:
        by_canon[canonicalize_name(tip)].append(tip)

    candidates = by_canon.get(canonicalize_name(foreground), [])
    if not candidates:
        raise RuntimeError(
            f"Could not map foreground branch '{foreground}' to any tree tip label. "
            f"Check naming consistency between aBSREL output and tree tip names."
        )
    if len(candidates) > 1:
        raise RuntimeError(
            f"Foreground '{foreground}' is ambiguous after normalization; matches: {', '.join(candidates)}"
        )
    return candidates[0]


def label_tree(tree: str, foreground: str) -> str:
    tree_label = resolve_foreground_name(tree, foreground)
    name = re.escape(tree_label)
    # PAML branch-site trees are most reliably parsed with foreground labels after branch length.
    labeled, n = re.subn(
        rf'(?:(?<=\()|(?<=,)){name}:([0-9eE.+-]+)',
        rf'{tree_label}:\1 #1',
        tree,
        count=1,
    )
    if n == 0:
        labeled = tree
    if labeled == tree:
        raise RuntimeError(
            f"Could not label foreground branch {foreground} (resolved tree label: {tree_label}) in tree"
        )
    return labeled


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument('--tree', required=True)
    p.add_argument('--foreground-list', required=True)
    p.add_argument('--outdir', required=True)
    p.add_argument('--manifest', required=True)
    a = p.parse_args()

    tree = Path(a.tree).read_text(encoding='utf-8').strip()
    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    foregrounds = [x.strip() for x in Path(a.foreground_list).read_text(encoding='utf-8').splitlines() if x.strip()]
    lines = []
    for fg in foregrounds:
        safe = fg.replace('/', '_')
        d = outdir / safe
        d.mkdir(parents=True, exist_ok=True)
        outfile = d / 'foreground.tree'
        outfile.write_text(label_tree(tree, fg) + '\n', encoding='utf-8')
        lines.append(f"{fg}\t{outfile}\n")
    Path(a.manifest).write_text(''.join(lines), encoding='utf-8')


if __name__ == '__main__':
    main()
