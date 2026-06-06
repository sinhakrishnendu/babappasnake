#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests


ALT_CTL = """      seqfile = {seqfile}
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
       model = 2
     NSsites = 2
       icode = 0
   fix_kappa = 0
       kappa = 2
   fix_omega = 0
       omega = 1.5
   fix_alpha = 1
       alpha = 0.
      Malpha = 0
       ncatG = 8
       getSE = 0
 RateAncestor = 0
   Small_Diff = .5e-6
    cleandata = 0
        method = 0
"""

NULL_CTL = ALT_CTL.replace("fix_omega = 0", "fix_omega = 1").replace("omega = 1.5", "omega = 1")


def parse_lnL(mlc_path: Path) -> float:
    for line in mlc_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        m = re.search(r"lnL\([^)]*\):\s*([-+0-9.eE]+)", line)
        if m:
            return float(m.group(1))
    raise RuntimeError(f"Could not parse lnL from {mlc_path}")


def run_codeml_and_get_lnl(codeml_exe: str, ctl_path: Path, workdir: Path, mlc_path: Path, foreground: str, model: str) -> float:
    if mlc_path.exists():
        mlc_path.unlink()

    res = subprocess.run(
        [codeml_exe, str(ctl_path.resolve())],
        cwd=workdir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        stdin=subprocess.DEVNULL,
    )

    try:
        lnl = parse_lnL(mlc_path)
    except Exception as exc:
        raise RuntimeError(
            f"codeml {model} failed for {foreground}\n"
            f"(exit code: {res.returncode})\n\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
        ) from exc

    # codeml can emit parse warnings to stderr and still produce valid results.
    if res.returncode != 0:
        print(
            f"[WARN] codeml {model} returned {res.returncode} for {foreground}, "
            f"but lnL was parsed successfully from {mlc_path}",
            file=sys.stderr,
        )
        if res.stderr.strip():
            print(f"[WARN] codeml {model} stderr for {foreground}: {res.stderr.strip()}", file=sys.stderr)

    return lnl


def run_one_foreground(fg: str, aln: str, tree_dir: Path, codeml_exe: str, codonfreq: int) -> tuple[str, float, float, float, float]:
    safe = fg.replace("/", "_")
    fg_dir = tree_dir / safe
    fg_dir.mkdir(parents=True, exist_ok=True)
    treefile = str((fg_dir / "foreground.tree").resolve())

    alt_dir = fg_dir / "alt"
    null_dir = fg_dir / "null"
    alt_dir.mkdir(exist_ok=True)
    null_dir.mkdir(exist_ok=True)

    alt_ctl = ALT_CTL.format(
        seqfile=aln,
        treefile=treefile,
        outfile=str((alt_dir / "mlc_alt.txt").resolve()),
        codonfreq=codonfreq,
    )
    null_ctl = NULL_CTL.format(
        seqfile=aln,
        treefile=treefile,
        outfile=str((null_dir / "mlc_null.txt").resolve()),
        codonfreq=codonfreq,
    )

    (alt_dir / "codeml.ctl").write_text(alt_ctl, encoding="utf-8")
    (null_dir / "codeml.ctl").write_text(null_ctl, encoding="utf-8")

    alt_lnl = run_codeml_and_get_lnl(
        codeml_exe,
        alt_dir / "codeml.ctl",
        alt_dir,
        alt_dir / "mlc_alt.txt",
        fg,
        "alt",
    )
    null_lnl = run_codeml_and_get_lnl(
        codeml_exe,
        null_dir / "codeml.ctl",
        null_dir,
        null_dir / "mlc_null.txt",
        fg,
        "null",
    )
    lrt = max(0.0, 2 * (alt_lnl - null_lnl))
    pval = chi2.sf(lrt, 1)
    return fg, alt_lnl, null_lnl, lrt, pval


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--alignment", required=True)
    p.add_argument("--tree-dir", required=True)
    p.add_argument("--foreground-list", required=True)
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--codeml", default="codeml")
    p.add_argument("--codonfreq", type=int, default=2)
    p.add_argument("--jobs", type=int, default=1)
    a = p.parse_args()

    aln = str(Path(a.alignment).resolve())
    tree_dir = Path(a.tree_dir)
    foregrounds = [x.strip() for x in Path(a.foreground_list).read_text(encoding="utf-8").splitlines() if x.strip()]
    rows_by_fg: dict[str, list[float | str]] = {}
    pvals_by_fg: dict[str, float] = {}

    jobs = max(1, int(a.jobs))
    if jobs == 1:
        for fg in foregrounds:
            fg_name, alt_lnl, null_lnl, lrt, pval = run_one_foreground(fg, aln, tree_dir, a.codeml, int(a.codonfreq))
            pvals_by_fg[fg_name] = pval
            rows_by_fg[fg_name] = [fg_name, alt_lnl, null_lnl, lrt, pval]
    else:
        with ThreadPoolExecutor(max_workers=jobs) as ex:
            fut_to_fg = {
                ex.submit(run_one_foreground, fg, aln, tree_dir, a.codeml, int(a.codonfreq)): fg
                for fg in foregrounds
            }
            for fut in as_completed(fut_to_fg):
                fg_name, alt_lnl, null_lnl, lrt, pval = fut.result()
                pvals_by_fg[fg_name] = pval
                rows_by_fg[fg_name] = [fg_name, alt_lnl, null_lnl, lrt, pval]

    rows = [rows_by_fg[fg] for fg in foregrounds if fg in rows_by_fg]
    pvals = [pvals_by_fg[fg] for fg in foregrounds if fg in pvals_by_fg]

    if rows:
        reject, qvals, _, _ = multipletests(pvals, alpha=0.05, method="fdr_bh")
    else:
        qvals = []
        reject = []

    with open(a.out_tsv, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["foreground_branch", "lnL_alt", "lnL_null", "LRT", "p_value", "q_value", "significant_BH_0.05"])
        for row, q, rej in zip(rows, qvals, reject):
            w.writerow(row + [q, str(bool(rej))])


if __name__ == "__main__":
    main()
