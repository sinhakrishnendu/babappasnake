#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

from Bio import Phylo


def branch_key(value):
    return re.sub(r"[^A-Za-z0-9]+", "_", (value or "")).strip("_")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--treefile", required=True)
    parser.add_argument("--branch-name", required=True)
    parser.add_argument("--output-tree", required=True)
    args = parser.parse_args()

    tree_path = Path(args.treefile)
    tree_text = tree_path.read_text().strip()
    tree = Phylo.read(args.treefile, "newick")
    matches = [
        clade
        for clade in tree.find_clades()
        if clade.name == args.branch_name or branch_key(clade.name) == args.branch_name
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one branch named {args.branch_name!r} in {args.treefile}, "
            f"found {len(matches)}."
        )
    target = matches[0]
    target_name = target.name
    if target.is_terminal():
        pattern = re.compile(rf"(?<![A-Za-z0-9_.-]){re.escape(target_name)}(?=[:),;])")
        rendered_tree, replacement_count = pattern.subn(f"{target_name}#1", tree_text, count=1)
    else:
        pattern = re.compile(rf"\){re.escape(target_name)}(?=[:),;])")
        rendered_tree, replacement_count = pattern.subn(")#1", tree_text, count=1)
    if replacement_count != 1:
        raise ValueError(
            f"Failed to render a codeml foreground tag for branch {args.branch_name!r} "
            f"in {args.treefile}."
        )
    Path(args.output_tree).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output_tree).write_text(rendered_tree + ("\n" if not rendered_tree.endswith("\n") else ""))


if __name__ == "__main__":
    main()
