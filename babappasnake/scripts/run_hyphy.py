#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path


def run(cmd: list[str]) -> None:
    res = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if res.returncode != 0:
        raise RuntimeError(
            "Command failed:\n"
            + " ".join(cmd)
            + f"\n\nExit code: {res.returncode}"
            + "\n\nSTDOUT:\n"
            + res.stdout
            + "\nSTDERR:\n"
            + res.stderr
        )


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--cds-aln", required=True)
    p.add_argument("--tree", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--hyphy", default="hyphy")
    a = p.parse_args()
    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    absrel_json = outdir / "absrel.json"
    meme_json = outdir / "meme.json"

    run([
        a.hyphy, "absrel",
        "--alignment", a.cds_aln,
        "--tree", a.tree,
        "--branches", "Leaves",
        "--output", str(absrel_json),
    ])
    run([
        a.hyphy, "meme",
        "--alignment", a.cds_aln,
        "--tree", a.tree,
        "--branches", "Leaves",
        "--output", str(meme_json),
    ])

    summary = {
        "hyphy_exe": a.hyphy,
        "absrel_json": str(absrel_json),
        "meme_json": str(meme_json),
    }
    (outdir / "hyphy_done.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
