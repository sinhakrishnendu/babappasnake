#!/usr/bin/env python3
from __future__ import annotations

import argparse
import difflib
import re
import shutil
from pathlib import Path

from Bio import Phylo


def normalize(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", text.lower())


def tokenize(text: str) -> list[str]:
    return [token for token in re.split(r"[^a-z0-9]+", text.lower()) if token]


def matches_outgroup(label: str, outgroup_query: str) -> bool:
    q = outgroup_query.strip().lower()
    l = label.strip().lower()
    if not q:
        return False
    if q in l:
        return True

    qn = normalize(q)
    ln = normalize(l)
    if qn and qn in ln:
        return True

    q_tokens = tokenize(q)
    l_tokens = tokenize(l)
    if q_tokens and all(any(qt in lt for lt in l_tokens) for qt in q_tokens):
        return True
    return False


def resolve_outgroup_matches(tree, outgroup_query: str):
    terminals = [terminal for terminal in tree.get_terminals() if terminal.name]
    matches = [terminal for terminal in terminals if matches_outgroup(terminal.name, outgroup_query)]
    if not matches:
        labels = [terminal.name for terminal in terminals]
        suggestions = difflib.get_close_matches(outgroup_query, labels, n=5, cutoff=0.3)
        hint = f" Suggestions: {suggestions}" if suggestions else ""
        raise RuntimeError(
            f"Could not match outgroup query '{outgroup_query}' to any tree tip labels.{hint}"
        )
    if len(matches) == len(terminals):
        raise RuntimeError(
            f"Outgroup query '{outgroup_query}' matched all tree tips; please provide a more specific label."
        )
    return matches


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--tree", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--outgroup", default="")
    a = p.parse_args()

    input_tree = Path(a.tree)
    output_tree = Path(a.output)
    output_tree.parent.mkdir(parents=True, exist_ok=True)

    outgroup_query = a.outgroup.strip()
    if not outgroup_query:
        shutil.copy2(input_tree, output_tree)
        return

    tree = Phylo.read(str(input_tree), "newick")
    matches = resolve_outgroup_matches(tree, outgroup_query)
    outgroup = matches[0] if len(matches) == 1 else tree.common_ancestor(matches)
    tree.root_with_outgroup(outgroup)
    Phylo.write(tree, str(output_tree), "newick")


if __name__ == "__main__":
    main()
