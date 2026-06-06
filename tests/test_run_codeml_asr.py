from __future__ import annotations

import json
import os
import stat
import subprocess
import sys
from pathlib import Path

import pytest


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def test_run_codeml_asr_does_not_accept_stale_rst_output(tmp_path):
    aln = tmp_path / "cds.aln.fasta"
    tree = tmp_path / "tree.nwk"
    outdir = tmp_path / "asr"
    fake_codeml = tmp_path / "fake_codeml.py"

    _write(aln, ">seq1\nATGAAAACCTAA\n>seq2\nATGAAAACCTAA\n")
    _write(tree, "(seq1:0.1,seq2:0.1);\n")
    _write(outdir / "rst", "stale rst\n")
    _write(outdir / "asr_done.json", json.dumps({"stale": True}))
    _write(
        fake_codeml,
        (
            "#!/usr/bin/env python3\n"
            "import pathlib, sys\n"
            "ctl = pathlib.Path(sys.argv[1]).read_text(encoding='utf-8')\n"
            "outfile = None\n"
            "for line in ctl.splitlines():\n"
            "    if line.strip().lower().startswith('outfile'):\n"
            "        outfile = line.split('=', 1)[1].strip()\n"
            "if outfile:\n"
            "    pathlib.Path(outfile).write_text('fresh mlc\\n', encoding='utf-8')\n"
            "sys.exit(0)\n"
        ),
    )
    os.chmod(fake_codeml, os.stat(fake_codeml).st_mode | stat.S_IEXEC)

    with pytest.raises(subprocess.CalledProcessError):
        subprocess.run(
            [
                sys.executable,
                "-m",
                "babappasnake.scripts.run_codeml_asr",
                "--alignment",
                str(aln),
                "--tree",
                str(tree),
                "--outdir",
                str(outdir),
                "--codeml",
                str(fake_codeml),
                "--codonfreq",
                "7",
            ],
            check=True,
        )

    assert (outdir / "mlc_asr.txt").exists()
    assert not (outdir / "rst").exists()
    assert not (outdir / "asr_done.json").exists()
