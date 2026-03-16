#!/usr/bin/env python3
import argparse
from pathlib import Path

from Bio import Phylo


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--treefile", required=True)
    parser.add_argument("--branch-name", required=True)
    parser.add_argument("--output-tree", required=True)
    args = parser.parse_args()

    tree = Phylo.read(args.treefile, "newick")
    matches = [clade for clade in tree.find_clades() if clade.name == args.branch_name]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one branch named {args.branch_name!r} in {args.treefile}, "
            f"found {len(matches)}."
        )
    target = matches[0]
    target.name = f"{target.name}#1"
    Path(args.output_tree).parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, args.output_tree, "newick")


if __name__ == "__main__":
    main()
