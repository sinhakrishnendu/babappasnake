#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import re


def canonicalize_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")


def extract_tip_names(tree: str) -> list[str]:
    return re.findall(r"(?:(?<=\()|(?<=,))([^():,]+)(?=:)", tree)


def resolve_foreground_name(tree: str, foreground: str) -> str:
    tip_names = extract_tip_names(tree)
    if foreground in tip_names:
        return foreground
    fg_canon = canonicalize_name(foreground)
    matches = [x for x in tip_names if canonicalize_name(x) == fg_canon]
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise RuntimeError(f"Could not map foreground '{foreground}' to any tree tip label")
    raise RuntimeError(f"Foreground '{foreground}' is ambiguous after normalization; matches: {', '.join(matches)}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--tree", required=True)
    p.add_argument("--foreground", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    tree = Path(a.tree).read_text(encoding="utf-8").strip()
    resolved = resolve_foreground_name(tree, a.foreground)
    name = re.escape(resolved)
    tree2, n = re.subn(
        rf'(?:(?<=\()|(?<=,)){name}:([0-9eE.+-]+)',
        rf'{resolved}:\1 #1',
        tree,
        count=1,
    )
    if n == 0:
        tree2 = tree
    if tree2 == tree:
        raise RuntimeError(f"Could not label resolved foreground '{resolved}' in tree")
    Path(a.out).write_text(tree2 + ("\n" if not tree2.endswith("\n") else ""), encoding="utf-8")


if __name__ == "__main__":
    main()
