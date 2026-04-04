from __future__ import annotations

import json
import os
import stat
import subprocess
import sys
from pathlib import Path


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def test_run_gard_screen_mode_none_writes_not_run_summary(tmp_path):
    cds_aln = tmp_path / "cds.aln.fasta"
    outdir = tmp_path / "recomb" / "m" / "t" / "gard"
    _write(cds_aln, ">a\nATGAAATAG\n>b\nATGAAATAG\n")

    subprocess.run(
        [
            sys.executable,
            "-m",
            "babappasnake.scripts.run_gard_screen",
            "--cds-aln",
            str(cds_aln),
            "--outdir",
            str(outdir),
            "--mode",
            "none",
        ],
        check=True,
    )

    summary = json.loads((outdir / "gard_summary.json").read_text(encoding="utf-8"))
    assert summary["status"] == "not_run"
    assert summary["completed"] is False
    assert summary["breakpoints_detected"] is False


def test_run_gard_screen_fallbacks_when_options_not_supported(tmp_path):
    cds_aln = tmp_path / "cds.aln.fasta"
    outdir = tmp_path / "recomb" / "m" / "t" / "gard"
    fake_hyphy = tmp_path / "hyphy_fake.py"
    _write(cds_aln, ">a\nATGAAATAG\n>b\nATGAAATAG\n")
    _write(
        fake_hyphy,
        (
            "#!/usr/bin/env python3\n"
            "import json, pathlib, sys\n"
            "args = sys.argv[1:]\n"
            "if '--mode' in args:\n"
            "    print('Unknown option: --mode', file=sys.stderr)\n"
            "    sys.exit(2)\n"
            "if '--output' not in args:\n"
            "    sys.exit(3)\n"
            "out = pathlib.Path(args[args.index('--output') + 1])\n"
            "out.parent.mkdir(parents=True, exist_ok=True)\n"
            "out.write_text(json.dumps({'breakpoints': [12, 44]}), encoding='utf-8')\n"
            "print('ok')\n"
            "sys.exit(0)\n"
        ),
    )
    os.chmod(fake_hyphy, os.stat(fake_hyphy).st_mode | stat.S_IEXEC)

    subprocess.run(
        [
            sys.executable,
            "-m",
            "babappasnake.scripts.run_gard_screen",
            "--cds-aln",
            str(cds_aln),
            "--outdir",
            str(outdir),
            "--hyphy",
            str(fake_hyphy),
            "--mode",
            "gard",
            "--gard-mode",
            "Faster",
            "--rate-classes",
            "3",
            "--fail-on-error",
        ],
        check=True,
    )

    summary = json.loads((outdir / "gard_summary.json").read_text(encoding="utf-8"))
    assert summary["status"] == "ok"
    assert summary["fallback_used"] is True
    assert summary["breakpoints_detected"] is True
    assert summary["n_breakpoints"] == 2

