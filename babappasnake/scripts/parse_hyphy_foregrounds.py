#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def extract_absrel_leaf_hits(data: dict, pcut: float) -> list[tuple[str, float]]:
    fits = data.get("branch attributes", {}).get("0", {})
    hits = []
    for branch, stats in fits.items():
        if not isinstance(stats, dict):
            continue
        is_leaf = stats.get("is leaf", None)
        p = stats.get("Corrected P-value", stats.get("p-value", stats.get("P-value")))
        if p is None:
            continue
        if (is_leaf is True or is_leaf is None) and float(p) <= pcut:
            hits.append((branch, float(p)))
    return sorted(hits, key=lambda x: x[1])


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--absrel-json", required=True)
    p.add_argument("--pcut", type=float, default=0.1)
    p.add_argument("--dynamic", action="store_true", help="Use dynamic thresholding: start at --dynamic-start and increase by --dynamic-step until at least one hit is found (or max reached).")
    p.add_argument("--dynamic-start", type=float, default=0.05)
    p.add_argument("--dynamic-step", type=float, default=0.01)
    p.add_argument("--dynamic-max", type=float, default=1.0)
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--out-list", required=True)
    p.add_argument("--out-meta", required=False)
    a = p.parse_args()

    data = json.loads(Path(a.absrel_json).read_text(encoding="utf-8"))
    selected_pcut = float(a.pcut)

    if a.dynamic:
        pcut = float(a.dynamic_start)
        max_p = float(a.dynamic_max)
        step = float(a.dynamic_step)
        if step <= 0:
            raise RuntimeError("--dynamic-step must be > 0")
        hits = []
        while pcut <= max_p + 1e-12:
            hits = extract_absrel_leaf_hits(data, pcut)
            selected_pcut = pcut
            if hits:
                break
            pcut = round(pcut + step, 10)
    else:
        hits = extract_absrel_leaf_hits(data, selected_pcut)

    with open(a.out_tsv, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["foreground_branch", "absrel_p"])
        w.writerows(hits)
    with open(a.out_list, "w", encoding="utf-8") as fh:
        for branch, _p in hits:
            fh.write(branch + "\n")

    if a.out_meta:
        meta = {
            "dynamic_mode": bool(a.dynamic),
            "selected_pcut": selected_pcut,
            "num_hits": len(hits),
            "dynamic_start": float(a.dynamic_start),
            "dynamic_step": float(a.dynamic_step),
            "dynamic_max": float(a.dynamic_max),
        }
        Path(a.out_meta).write_text(json.dumps(meta, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
