#!/usr/bin/env python3
from __future__ import annotations

import argparse
import difflib
import re
import shutil
import sys
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


def apply_outgroup_or_copy(input_tree: Path, output_tree: Path, outgroup_query: str) -> dict[str, str]:
    output_tree.parent.mkdir(parents=True, exist_ok=True)
    normalized_query = outgroup_query.strip()
    if not normalized_query:
        shutil.copy2(input_tree, output_tree)
        return {"status": "copied_unrooted", "reason": "no_outgroup_supplied"}

    tree = Phylo.read(str(input_tree), "newick")
    try:
        matches = resolve_outgroup_matches(tree, normalized_query)
        outgroup = matches[0] if len(matches) == 1 else tree.common_ancestor(matches)
        tree.root_with_outgroup(outgroup)
        Phylo.write(tree, str(output_tree), "newick")
        return {"status": "rooted", "reason": f"matched_{len(matches)}_tip(s)"}
    except Exception as exc:
        shutil.copy2(input_tree, output_tree)
        print(
            f"[WARN] Could not apply outgroup '{normalized_query}'; using unrooted tree downstream. {exc}",
            file=sys.stderr,
        )
        return {"status": "copied_unrooted", "reason": "outgroup_not_applied"}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--tree", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--outgroup", nargs="?", const="", default="")
    a = p.parse_args()

    input_tree = Path(a.tree)
    output_tree = Path(a.output)
    apply_outgroup_or_copy(input_tree, output_tree, str(a.outgroup))


if __name__ == "__main__":
    main()
