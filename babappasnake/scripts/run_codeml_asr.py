#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path


ASR_CTL = """      seqfile = {seqfile}
     treefile = {treefile}
      outfile = {outfile}
        noisy = 3
      verbose = 1
      runmode = 0
      seqtype = 1
    CodonFreq = {codonfreq}
        clock = 0
       aaDist = 0
           aaRatefile = wag.dat
       model = 0
     NSsites = 0
       icode = 0
   fix_kappa = 0
       kappa = 2
   fix_omega = 0
       omega = 1
   fix_alpha = 1
       alpha = 0.
      Malpha = 0
       ncatG = 8
       getSE = 0
 RateAncestor = 1
   Small_Diff = .5e-6
    cleandata = 0
        method = 0
"""


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--alignment", required=True)
    p.add_argument("--tree", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--codeml", default="codeml")
    p.add_argument("--codonfreq", type=int, default=2)
    a = p.parse_args()

    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    mlc = outdir / "mlc_asr.txt"
    rst = outdir / "rst"
    ctl = outdir / "codeml_asr.ctl"
    ctl.write_text(
        ASR_CTL.format(
            seqfile=str(Path(a.alignment).resolve()),
            treefile=str(Path(a.tree).resolve()),
            outfile=str(mlc.resolve()),
            codonfreq=int(a.codonfreq),
        ),
        encoding="utf-8",
    )

    if mlc.exists():
        mlc.unlink()
    if rst.exists():
        rst.unlink()
    done_json = outdir / "asr_done.json"
    if done_json.exists():
        done_json.unlink()

    res = subprocess.run(
        [a.codeml, str(ctl.resolve())],
        cwd=outdir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        stdin=subprocess.DEVNULL,
    )

    # codeml may return non-zero while still producing valid outputs.
    if not mlc.exists() or not rst.exists():
        raise RuntimeError(
            f"codeml ASR failed (exit code: {res.returncode})\n\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
        )

    summary = {
        "codeml": a.codeml,
        "mlc_asr": str(mlc.resolve()),
        "rst": str(rst.resolve()),
    }
    done_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
