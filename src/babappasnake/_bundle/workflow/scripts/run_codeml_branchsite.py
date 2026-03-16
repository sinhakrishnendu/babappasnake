#!/usr/bin/env python3
import argparse
import shutil
import subprocess
from pathlib import Path


def write_ctl(path, seqfile, treefile, outfile, fix_omega, omega, args):
    path.write_text(
        "\n".join(
            [
                f"seqfile = {seqfile}",
                f"treefile = {treefile}",
                f"outfile = {outfile}",
                f"noisy = {args.noisy}",
                f"verbose = {args.verbose}",
                "runmode = 0",
                "seqtype = 1",
                f"CodonFreq = {args.codon_frequency}",
                "clock = 0",
                "aaDist = 0",
                "model = 2",
                "NSsites = 2",
                f"cleandata = {args.cleandata}",
                "fix_kappa = 0",
                f"kappa = {args.kappa_initial}",
                "fix_alpha = 1",
                "alpha = 0",
                "Malpha = 0",
                "ncatG = 8",
                "getSE = 0",
                "RateAncestor = 0",
                "method = 0",
                f"fix_omega = {fix_omega}",
                f"omega = {omega}",
            ]
        )
        + "\n"
    )


def run_model(codeml, ctl_path, workdir, stdout_log):
    with open(stdout_log, "w") as handle:
        result = subprocess.run([codeml, ctl_path.name], cwd=workdir, stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        raise RuntimeError(f"codeml failed for {ctl_path} with exit code {result.returncode}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--codeml", required=True)
    parser.add_argument("--alignment-fasta", required=True)
    parser.add_argument("--foreground-tree", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--branch-name", required=True)
    parser.add_argument("--codon-frequency", required=True, type=int)
    parser.add_argument("--cleandata", required=True, type=int)
    parser.add_argument("--kappa-initial", required=True, type=float)
    parser.add_argument("--omega-initial", required=True, type=float)
    parser.add_argument("--noisy", required=True, type=int)
    parser.add_argument("--verbose", required=True, type=int)
    parser.add_argument("--done-out", required=True)
    args = parser.parse_args()

    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    alignment_copy = workdir / "aligned_cds.fasta"
    tree_copy = workdir / "foreground.treefile"
    shutil.copyfile(args.alignment_fasta, alignment_copy)
    shutil.copyfile(args.foreground_tree, tree_copy)

    alt_dir = workdir / "branch_site_alt"
    null_dir = workdir / "branch_site_null"
    alt_dir.mkdir(parents=True, exist_ok=True)
    null_dir.mkdir(parents=True, exist_ok=True)

    alt_ctl = alt_dir / "codeml.ctl"
    null_ctl = null_dir / "codeml.ctl"
    write_ctl(alt_ctl, "../aligned_cds.fasta", "../foreground.treefile", "output.txt", 0, args.omega_initial, args)
    write_ctl(null_ctl, "../aligned_cds.fasta", "../foreground.treefile", "output.txt", 1, 1.0, args)

    run_model(args.codeml, alt_ctl, alt_dir, alt_dir / "codeml.stdout.log")
    run_model(args.codeml, null_ctl, null_dir, null_dir / "codeml.stdout.log")

    if not (alt_dir / "output.txt").exists() or not (null_dir / "output.txt").exists():
        raise FileNotFoundError(
            f"codeml did not produce both branch-site outputs for branch {args.branch_name}"
        )

    Path(args.done_out).write_text("completed\n")


if __name__ == "__main__":
    main()
