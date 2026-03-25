from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def test_run_hyphy_writes_fallback_json_when_hyphy_fails(tmp_path):
    cds_aln = tmp_path / "cds.aln.fasta"
    tree = tmp_path / "tree.nwk"
    outdir = tmp_path / "hyphy"

    _write(cds_aln, ">seq1\nATGAAAACCTAA\n>seq2\nATGAAAACCTAA\n")
    _write(tree, "(seq1:0.1,seq2:0.1);\n")

    subprocess.run(
        [
            sys.executable,
            "-m",
            "babappasnake.scripts.run_hyphy",
            "--cds-aln",
            str(cds_aln),
            "--tree",
            str(tree),
            "--outdir",
            str(outdir),
            "--hyphy",
            "/usr/bin/false",
            "--absrel-branches",
            "Leaves",
            "--meme-branches",
            "Leaves",
        ],
        check=True,
    )

    absrel = json.loads((outdir / "absrel.json").read_text(encoding="utf-8"))
    meme = json.loads((outdir / "meme.json").read_text(encoding="utf-8"))
    done = json.loads((outdir / "hyphy_done.json").read_text(encoding="utf-8"))

    assert absrel.get("status") == "failed"
    assert meme.get("status") == "failed"
    assert done.get("absrel_status") == "failed"
    assert done.get("meme_status") == "failed"
    assert "branch attributes" in absrel
