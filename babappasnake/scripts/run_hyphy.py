#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path


def run_command(cmd: list[str]) -> tuple[int, str, str]:
    try:
        res = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return int(res.returncode), res.stdout, res.stderr
    except Exception as exc:
        return 127, "", f"{type(exc).__name__}: {exc}"


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--cds-aln", required=True)
    p.add_argument("--tree", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--hyphy", default="hyphy")
    p.add_argument("--absrel-branches", default="Leaves")
    p.add_argument("--meme-branches", default="Leaves")
    p.add_argument("--fail-on-error", action="store_true", help="Fail immediately if a HyPhy command fails.")
    a = p.parse_args()
    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    absrel_json = outdir / "absrel.json"
    meme_json = outdir / "meme.json"

    absrel_cmd = [
        a.hyphy, "absrel",
        "--alignment", a.cds_aln,
        "--tree", a.tree,
        "--branches", a.absrel_branches,
        "--output", str(absrel_json),
    ]
    meme_cmd = [
        a.hyphy, "meme",
        "--alignment", a.cds_aln,
        "--tree", a.tree,
        "--branches", a.meme_branches,
        "--output", str(meme_json),
    ]

    absrel_rc, absrel_stdout, absrel_stderr = run_command(absrel_cmd)
    meme_rc, meme_stdout, meme_stderr = run_command(meme_cmd)

    if absrel_rc != 0:
        write_json(
            absrel_json,
            {
                "status": "failed",
                "analysis": "absrel",
                "command": absrel_cmd,
                "return_code": absrel_rc,
                "stdout": absrel_stdout,
                "stderr": absrel_stderr,
                "branch attributes": {"0": {}},
            },
        )
    elif not absrel_json.exists() or absrel_json.stat().st_size == 0:
        write_json(
            absrel_json,
            {
                "status": "failed",
                "analysis": "absrel",
                "command": absrel_cmd,
                "return_code": 125,
                "stdout": absrel_stdout,
                "stderr": absrel_stderr,
                "error": "HyPhy absrel finished but output JSON is missing or empty.",
                "branch attributes": {"0": {}},
            },
        )
        absrel_rc = 125

    if meme_rc != 0:
        write_json(
            meme_json,
            {
                "status": "failed",
                "analysis": "meme",
                "command": meme_cmd,
                "return_code": meme_rc,
                "stdout": meme_stdout,
                "stderr": meme_stderr,
            },
        )
    elif not meme_json.exists() or meme_json.stat().st_size == 0:
        write_json(
            meme_json,
            {
                "status": "failed",
                "analysis": "meme",
                "command": meme_cmd,
                "return_code": 125,
                "stdout": meme_stdout,
                "stderr": meme_stderr,
                "error": "HyPhy meme finished but output JSON is missing or empty.",
            },
        )
        meme_rc = 125

    if a.fail_on_error and (absrel_rc != 0 or meme_rc != 0):
        raise RuntimeError(
            f"HyPhy failed (absrel_rc={absrel_rc}, meme_rc={meme_rc}). "
            f"See {absrel_json} and {meme_json} for details."
        )

    summary = {
        "hyphy_exe": a.hyphy,
        "absrel_json": str(absrel_json),
        "meme_json": str(meme_json),
        "absrel_status": "ok" if absrel_rc == 0 else "failed",
        "meme_status": "ok" if meme_rc == 0 else "failed",
        "absrel_return_code": absrel_rc,
        "meme_return_code": meme_rc,
    }
    write_json(outdir / "hyphy_done.json", summary)


if __name__ == "__main__":
    main()
