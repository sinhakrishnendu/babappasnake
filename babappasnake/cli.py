from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path
import subprocess
import yaml
import importlib.resources as ir

from babappasnake.utils import format_missing_tools, resolve_tools


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="babappasnake",
        description="A simple Snakemake pipeline for episodic selection analysis.",
    )
    p.add_argument("--prot", required=True, help="Path to proteomes folder")
    p.add_argument("--query", required=True, help="Path to protein query FASTA")
    p.add_argument("--cds", default=None, help="Optional CDS FASTA for the orthogroup")
    p.add_argument("--outdir", default="babappasnake_run", help="Output directory")
    p.add_argument("--coverage", type=float, default=0.70, help="RBH minimum reciprocal coverage [default: 0.70]")
    p.add_argument("--threads", type=int, default=4, help="Threads for external tools [default: 4]")
    p.add_argument("--absrel-p", type=float, default=0.1, help="Leaf-branch significance threshold for aBSREL [default: 0.1]")
    p.add_argument("--meme-p", type=float, default=0.1, help="MEME site threshold retained in summary [default: 0.1]")
    p.add_argument("--use-clipkit", choices=["yes", "no"], default="yes", help="Trim alignments with ClipKIT [default: yes]")
    p.add_argument("--snake-args", default="", help="Extra raw arguments forwarded to snakemake")
    return p.parse_args()


def stage_inputs(args: argparse.Namespace, outdir: Path) -> None:
    staged = outdir / "inputs"
    staged.mkdir(parents=True, exist_ok=True)
    shutil.copy2(args.query, staged / "query.fasta")
    prot_dst = staged / "proteomes"
    prot_dst.mkdir(exist_ok=True)
    for item in sorted(Path(args.prot).iterdir()):
        if item.is_file():
            shutil.copy2(item, prot_dst / item.name)
    if args.cds:
        cds_dst_dir = outdir / "user_supplied"
        cds_dst_dir.mkdir(exist_ok=True)
        shutil.copy2(args.cds, cds_dst_dir / "orthogroup_cds.fasta")


def write_config(args: argparse.Namespace, outdir: Path, executables: dict[str, str]) -> Path:
    cfg = {
        "outdir": str(outdir.resolve()),
        "query_fasta": str((outdir / "inputs" / "query.fasta").resolve()),
        "proteomes_dir": str((outdir / "inputs" / "proteomes").resolve()),
        "user_cds": str((outdir / "user_supplied" / "orthogroup_cds.fasta").resolve()),
        "coverage": float(args.coverage),
        "threads": int(args.threads),
        "absrel_p": float(args.absrel_p),
        "meme_p": float(args.meme_p),
        "use_clipkit": args.use_clipkit == "yes",
        "executables": executables,
    }
    config_path = outdir / "config.yaml"
    with open(config_path, "w", encoding="utf-8") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)
    return config_path


def main() -> None:
    args = parse_args()
    resolved, missing = resolve_tools()
    if missing:
        print(format_missing_tools(missing), file=sys.stderr)
        raise SystemExit(2)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stage_inputs(args, outdir)
    config_path = write_config(args, outdir, resolved)

    snakefile = ir.files("babappasnake").joinpath("workflow", "Snakefile")
    cmd = [
        sys.executable,
        "-m",
        "snakemake",
        "--snakefile",
        str(snakefile),
        "--configfile",
        str(config_path),
        "--cores",
        str(args.threads),
        "--printshellcmds",
        "--rerun-incomplete",
    ]
    if args.snake_args.strip():
        cmd.extend(args.snake_args.strip().split())

    result = subprocess.run(cmd, check=False)
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
