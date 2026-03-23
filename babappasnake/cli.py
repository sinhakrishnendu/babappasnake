from __future__ import annotations

import argparse
import importlib.resources as ir
import json
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import yaml

from babappasnake.utils import format_missing_tools, resolve_tools


IQTREE_BOOTSTRAP_CHOICES = (1000, 5000, 10000)
CLIPKIT_MODE_CHOICES = (
    "kpi",
    "kpic",
    "gappy",
    "smart-gap",
    "kpi-gappy",
    "kpic-gappy",
    "kpi-smart-gap",
    "kpic-smart-gap",
    "strict",
    "relaxed",
)


@dataclass(frozen=True)
class StepSpec:
    rule: str
    description: str
    outputs: tuple[str, ...]


def is_tty_interactive() -> bool:
    return sys.stdin.isatty() and sys.stdout.isatty()


def prompt_text(label: str, default: str = "", required: bool = False) -> str:
    while True:
        shown_default = f" [{default}]" if default else ""
        value = input(f"{label}{shown_default}: ").strip()
        if value:
            return value
        if default:
            return default
        if not required:
            return ""
        print("This field is required.")


def prompt_yes_no(label: str, default: bool = True) -> bool:
    suffix = " [Y/n]" if default else " [y/N]"
    while True:
        value = input(f"{label}{suffix}: ").strip().lower()
        if not value:
            return default
        if value in {"y", "yes"}:
            return True
        if value in {"n", "no"}:
            return False
        print("Please answer y or n.")


def prompt_float(label: str, default: float) -> float:
    while True:
        value = input(f"{label} [{default}]: ").strip()
        if not value:
            return float(default)
        try:
            return float(value)
        except ValueError:
            print("Please enter a numeric value.")


def prompt_int(label: str, default: int) -> int:
    while True:
        value = input(f"{label} [{default}]: ").strip()
        if not value:
            return int(default)
        try:
            return int(value)
        except ValueError:
            print("Please enter an integer value.")


def prompt_choice(label: str, choices: tuple[str, ...], default: str) -> str:
    print(label)
    for idx, option in enumerate(choices, start=1):
        mark = " (default)" if option == default else ""
        print(f"  {idx}. {option}{mark}")
    while True:
        value = input("Select option number or press Enter for default: ").strip()
        if not value:
            return default
        if value.isdigit():
            idx = int(value)
            if 1 <= idx <= len(choices):
                return choices[idx - 1]
        if value in choices:
            return value
        print("Invalid selection.")


def prompt_step_action(can_skip: bool) -> str:
    options = "run / stop" if not can_skip else "run / skip / stop"
    while True:
        value = input(f"Choose action ({options}) [run]: ").strip().lower()
        if not value or value == "run":
            return "run"
        if can_skip and value == "skip":
            return "skip"
        if value in {"stop", "quit", "exit"}:
            return "stop"
        print("Invalid selection.")


def prompt_bootstrap(default: int) -> int:
    print("IQ-TREE bootstrap replicates")
    for idx, option in enumerate(IQTREE_BOOTSTRAP_CHOICES, start=1):
        mark = " (default)" if option == default else ""
        print(f"  {idx}. {option}{mark}")
    print("  c. custom integer")
    while True:
        value = input("Select option (1/2/3/c) or press Enter for default: ").strip().lower()
        if not value:
            return default
        if value in {"1", "2", "3"}:
            return IQTREE_BOOTSTRAP_CHOICES[int(value) - 1]
        if value == "c":
            return prompt_int("Enter custom bootstrap replicates", default)
        if value.isdigit():
            return int(value)
        print("Invalid selection.")


def maybe_prompt_interactive(args: argparse.Namespace) -> argparse.Namespace:
    if args.interactive == "no":
        return args
    if not is_tty_interactive():
        return args

    print("\nInteractive babappasnake configuration\n")
    args.prot = prompt_text("Proteomes directory (--prot)", args.prot or "", required=True)
    args.query = prompt_text("Query FASTA (--query)", args.query or "", required=True)
    args.cds = prompt_text("CDS FASTA (--cds, optional)", args.cds or "", required=False) or None
    args.outgroup = prompt_text("Outgroup query (--outgroup, optional)", args.outgroup or "", required=False)
    args.outdir = prompt_text("Output directory (--outdir)", args.outdir, required=True)
    args.coverage = prompt_float("RBH reciprocal coverage (--coverage)", float(args.coverage))
    args.threads = prompt_int("Threads (--threads)", int(args.threads))
    args.use_clipkit = "yes" if prompt_yes_no("Use ClipKIT?", args.use_clipkit == "yes") else "no"
    if args.use_clipkit == "yes":
        args.clipkit_mode_protein = prompt_choice(
            "ClipKIT protein mode (--clipkit-mode-protein)",
            CLIPKIT_MODE_CHOICES,
            args.clipkit_mode_protein,
        )
        args.clipkit_mode_codon = prompt_choice(
            "ClipKIT codon mode (--clipkit-mode-codon)",
            CLIPKIT_MODE_CHOICES,
            args.clipkit_mode_codon,
        )
    args.iqtree_bootstrap = prompt_bootstrap(int(args.iqtree_bootstrap))
    args.iqtree_bnni = "yes" if prompt_yes_no("Use IQ-TREE -bnni?", args.iqtree_bnni == "yes") else "no"
    args.iqtree_model = prompt_text("IQ-TREE model string (--iqtree-model)", args.iqtree_model, required=True)
    args.absrel_p = prompt_float("aBSREL compatibility threshold (--absrel-p)", float(args.absrel_p))
    args.absrel_dynamic_start = prompt_float(
        "Dynamic aBSREL start p (--absrel-dynamic-start)",
        float(args.absrel_dynamic_start),
    )
    args.absrel_dynamic_step = prompt_float(
        "Dynamic aBSREL increment (--absrel-dynamic-step)",
        float(args.absrel_dynamic_step),
    )
    args.absrel_dynamic_max = prompt_float(
        "Dynamic aBSREL max p (--absrel-dynamic-max)",
        float(args.absrel_dynamic_max),
    )
    args.meme_p = prompt_float("MEME report threshold (--meme-p)", float(args.meme_p))
    args.snake_args = prompt_text("Extra snakemake args (--snake-args, optional)", args.snake_args or "", required=False)
    args.guided = "yes" if prompt_yes_no("Run in guided step-by-step mode?", args.guided == "yes") else "no"

    print("\nConfiguration summary:")
    print(f"- prot: {args.prot}")
    print(f"- query: {args.query}")
    print(f"- cds: {args.cds or '(none)'}")
    print(f"- outgroup: {args.outgroup or '(none)'}")
    print(f"- outdir: {args.outdir}")
    print(f"- threads: {args.threads}")
    print(f"- coverage: {args.coverage}")
    print(f"- use_clipkit: {args.use_clipkit}")
    if args.use_clipkit == "yes":
        print(f"- clipkit_mode_protein: {args.clipkit_mode_protein}")
        print(f"- clipkit_mode_codon: {args.clipkit_mode_codon}")
    print(f"- iqtree_bootstrap: {args.iqtree_bootstrap}")
    print(f"- iqtree_bnni: {args.iqtree_bnni}")
    print(f"- iqtree_model: {args.iqtree_model}")
    print(f"- guided: {args.guided}")
    if not prompt_yes_no("Proceed with this configuration?", True):
        raise SystemExit(1)
    print("")
    return args


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="babappasnake",
        description="A simple Snakemake pipeline for episodic selection analysis.",
    )
    p.add_argument("--prot", default=None, help="Path to proteomes folder")
    p.add_argument("--query", default=None, help="Path to protein query FASTA")
    p.add_argument("--cds", default=None, help="Optional CDS FASTA for the orthogroup")
    p.add_argument("--outdir", default="babappasnake_run", help="Output directory")
    p.add_argument("--coverage", type=float, default=0.70, help="RBH minimum reciprocal coverage [default: 0.70]")
    p.add_argument("--threads", type=int, default=4, help="Threads for external tools [default: 4]")
    p.add_argument("--outgroup", default="", help="Outgroup label query used to root the IQ-TREE output (case-insensitive substring match)")
    p.add_argument("--iqtree-bootstrap", type=int, default=1000, help="IQ-TREE ultrafast bootstrap replicates [default: 1000]")
    p.add_argument("--iqtree-bnni", choices=["yes", "no"], default="no", help="Enable IQ-TREE -bnni [default: no]")
    p.add_argument("--iqtree-model", default="MFP", help="IQ-TREE model string [default: MFP]")
    p.add_argument("--clipkit-mode-protein", default="kpic-smart-gap", help="ClipKIT mode for protein trimming")
    p.add_argument("--clipkit-mode-codon", default="kpic-smart-gap", help="ClipKIT mode for codon trimming")
    p.add_argument("--absrel-p", type=float, default=0.1, help="Leaf-branch significance threshold for aBSREL [default: 0.1]")
    p.add_argument("--absrel-dynamic-start", type=float, default=0.05, help="Dynamic aBSREL start threshold [default: 0.05]")
    p.add_argument("--absrel-dynamic-step", type=float, default=0.01, help="Dynamic aBSREL step size [default: 0.01]")
    p.add_argument("--absrel-dynamic-max", type=float, default=1.0, help="Dynamic aBSREL max threshold [default: 1.0]")
    p.add_argument("--meme-p", type=float, default=0.1, help="MEME site threshold retained in summary [default: 0.1]")
    p.add_argument("--use-clipkit", choices=["yes", "no"], default="yes", help="Trim alignments with ClipKIT [default: yes]")
    p.add_argument("--snake-args", default="", help="Extra raw arguments forwarded to snakemake")
    p.add_argument("--interactive", choices=["yes", "no"], default="yes", help="Prompt for pipeline settings interactively [default: yes]")
    p.add_argument("--guided", choices=["yes", "no"], default="yes", help="Run step-by-step guided execution [default: yes]")
    return p.parse_args()


def validate_inputs(args: argparse.Namespace) -> None:
    missing = []
    if not args.prot:
        missing.append("--prot")
    if not args.query:
        missing.append("--query")
    if missing:
        raise SystemExit(f"Missing required arguments: {', '.join(missing)}")

    prot_path = Path(args.prot)
    if not prot_path.is_dir():
        raise SystemExit(f"--prot must be an existing directory: {prot_path}")
    query_path = Path(args.query)
    if not query_path.is_file():
        raise SystemExit(f"--query must be an existing FASTA file: {query_path}")
    if args.cds:
        cds_path = Path(args.cds)
        if not cds_path.is_file():
            raise SystemExit(f"--cds must point to an existing FASTA file: {cds_path}")


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
        "outgroup_query": str(args.outgroup).strip(),
        "iqtree_bootstrap": int(args.iqtree_bootstrap),
        "iqtree_bnni": args.iqtree_bnni == "yes",
        "iqtree_model": str(args.iqtree_model).strip(),
        "use_clipkit": args.use_clipkit == "yes",
        "clipkit_mode_protein": str(args.clipkit_mode_protein).strip(),
        "clipkit_mode_codon": str(args.clipkit_mode_codon).strip(),
        "absrel_p": float(args.absrel_p),
        "absrel_dynamic_start": float(args.absrel_dynamic_start),
        "absrel_dynamic_step": float(args.absrel_dynamic_step),
        "absrel_dynamic_max": float(args.absrel_dynamic_max),
        "meme_p": float(args.meme_p),
        "executables": executables,
    }
    config_path = outdir / "config.yaml"
    with open(config_path, "w", encoding="utf-8") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)
    return config_path


def build_snakemake_cmd(config_path: Path, cores: int, target: str, snake_args: str, snakefile: Path) -> list[str]:
    cmd = [
        sys.executable,
        "-m",
        "snakemake",
        "--snakefile",
        str(snakefile),
        "--configfile",
        str(config_path),
        "--cores",
        str(cores),
        "--printshellcmds",
        "--rerun-incomplete",
        target,
    ]
    if snake_args.strip():
        try:
            cmd.extend(shlex.split(snake_args))
        except ValueError as exc:
            raise SystemExit(f"Invalid --snake-args value: {exc}") from exc
    return cmd


def run_snakemake_target(config_path: Path, cores: int, target: str, snake_args: str, snakefile: Path) -> int:
    cmd = build_snakemake_cmd(config_path, cores, target, snake_args, snakefile)
    result = subprocess.run(cmd, check=False)
    return result.returncode


def build_step_plan(have_cds: bool, use_clipkit: bool) -> list[StepSpec]:
    if not have_cds:
        return [
            StepSpec(
                "rbh_orthogroup",
                "Build orthogroup proteins and RBH mapping.",
                (
                    "orthogroup/orthogroup_proteins.fasta",
                    "orthogroup/orthogroup_headers.txt",
                    "orthogroup/rbh_summary.tsv",
                ),
            ),
            StepSpec(
                "waiting_note",
                "Write waiting note for later CDS submission.",
                ("orthogroup/WAITING_FOR_CDS.txt",),
            ),
        ]

    steps: list[StepSpec] = [
        StepSpec(
            "rbh_orthogroup",
            "Build orthogroup proteins and RBH mapping.",
            (
                "orthogroup/orthogroup_proteins.fasta",
                "orthogroup/orthogroup_headers.txt",
                "orthogroup/rbh_summary.tsv",
            ),
        ),
        StepSpec(
            "map_cds",
            "Map user CDS to orthogroup proteins with CDS quality filtering.",
            (
                "mapped_cds/mapped_orthogroup_cds.fasta",
                "mapped_cds/mapped_orthogroup_proteins.fasta",
                "mapped_cds/cds_protein_mapping.tsv",
            ),
        ),
        StepSpec(
            "align_proteins",
            "Generate protein alignment.",
            ("alignments/orthogroup_proteins.protein.aln.fasta",),
        ),
        StepSpec(
            "align_cds",
            "Generate codon alignment.",
            (
                "alignments/mapped_orthogroup_cds.protein.aln.fasta",
                "alignments/mapped_orthogroup_cds.codon.aln.fasta",
            ),
        ),
    ]
    if use_clipkit:
        steps.extend(
            [
                StepSpec(
                    "trim_protein_alignment",
                    "Trim protein alignment with ClipKIT.",
                    ("trimmed/orthogroup_proteins.clipkit.fasta",),
                ),
                StepSpec(
                    "trim_codon_alignment",
                    "Trim codon alignment with ClipKIT.",
                    ("trimmed/mapped_orthogroup_cds.clipkit.fasta",),
                ),
                StepSpec(
                    "strip_terminal_stop_codon",
                    "Remove terminal stop codons from trimmed codon alignment.",
                    ("trimmed/mapped_orthogroup_cds.clipkit.nostop.fasta",),
                ),
            ]
        )
    steps.extend(
        [
            StepSpec("iqtree_ml", "Infer ML phylogeny with IQ-TREE.", ("tree/orthogroup.treefile",)),
            StepSpec(
                "root_iqtree_outgroup",
                "Root the tree using the outgroup query.",
                ("tree/orthogroup.rooted.treefile",),
            ),
            StepSpec(
                "hyphy_exploratory",
                "Run HyPhy aBSREL + MEME exploratory tests.",
                ("hyphy/hyphy_done.json", "hyphy/absrel.json", "hyphy/meme.json"),
            ),
            StepSpec(
                "parse_foregrounds",
                "Select significant foreground branches.",
                (
                    "hyphy/significant_foregrounds.tsv",
                    "hyphy/significant_foregrounds.txt",
                    "hyphy/foreground_threshold.json",
                ),
            ),
            StepSpec(
                "prepare_foreground_trees",
                "Prepare branch-labeled trees for branch-site codeml.",
                ("branchsite/foreground_trees.tsv",),
            ),
            StepSpec(
                "branchsite_batch",
                "Run branch-site codeml across selected foregrounds.",
                ("branchsite/branchsite_results.tsv",),
            ),
            StepSpec(
                "codeml_asr",
                "Run codeml ancestral sequence reconstruction.",
                ("asr/asr_done.json", "asr/mlc_asr.txt", "asr/rst"),
            ),
            StepSpec(
                "final_summary",
                "Write final episodic selection summary.",
                ("summary/episodic_selection_summary.txt",),
            ),
        ]
    )
    return steps


def print_step_outputs(outdir: Path, outputs: tuple[str, ...]) -> None:
    def _print_preview(path: Path) -> None:
        suffix = path.suffix.lower()
        preview_extensions = {".txt", ".tsv", ".json", ".treefile", ".fasta", ".fa", ".aln"}
        if suffix not in preview_extensions:
            return
        try:
            if suffix == ".json":
                with open(path, "r", encoding="utf-8") as fh:
                    payload = json.load(fh)
                compact = json.dumps(payload, indent=2)
                lines = compact.splitlines()[:10]
            else:
                with open(path, "r", encoding="utf-8") as fh:
                    lines = [line.rstrip("\n") for _, line in zip(range(10), fh)]
            if not lines:
                return
            print("    preview:")
            for line in lines:
                print(f"      {line}")
        except Exception:
            return

    print("Step outputs:")
    for rel in outputs:
        path = outdir / rel
        if path.exists():
            if path.is_file():
                print(f"  - {rel} (file, {path.stat().st_size} bytes)")
                _print_preview(path)
            else:
                print(f"  - {rel} (directory)")
        else:
            print(f"  - {rel} (missing)")


def run_guided_pipeline(
    config_path: Path,
    outdir: Path,
    cores: int,
    use_clipkit: bool,
    snake_args: str,
    snakefile: Path,
) -> int:
    user_cds = outdir / "user_supplied" / "orthogroup_cds.fasta"
    have_cds = user_cds.exists() and user_cds.stat().st_size > 0
    steps = build_step_plan(have_cds, use_clipkit)
    print(f"Guided mode enabled. Planned steps: {len(steps)}")

    for idx, step in enumerate(steps, start=1):
        print(f"\nStep {idx}/{len(steps)}: {step.rule}")
        print(step.description)
        existing = tuple((outdir / rel).exists() for rel in step.outputs)
        can_skip = all(existing)
        if can_skip:
            print("All expected outputs for this step already exist.")
        action = prompt_step_action(can_skip=can_skip)
        if action == "stop":
            print("Guided run stopped by user.")
            return 130
        if action == "skip":
            print("Step skipped (existing outputs retained).")
            print_step_outputs(outdir, step.outputs)
            continue
        code = run_snakemake_target(config_path, cores, step.rule, snake_args, snakefile)
        if code != 0:
            return code
        print_step_outputs(outdir, step.outputs)
    return 0


def main() -> None:
    args = parse_args()
    args = maybe_prompt_interactive(args)
    validate_inputs(args)

    resolved, missing = resolve_tools()
    if missing:
        print(format_missing_tools(missing), file=sys.stderr)
        raise SystemExit(2)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stage_inputs(args, outdir)
    config_path = write_config(args, outdir, resolved)
    snakefile = ir.files("babappasnake").joinpath("workflow", "Snakefile")

    if args.guided == "yes" and is_tty_interactive():
        rc = run_guided_pipeline(
            config_path=config_path,
            outdir=outdir,
            cores=int(args.threads),
            use_clipkit=args.use_clipkit == "yes",
            snake_args=args.snake_args,
            snakefile=snakefile,
        )
    else:
        rc = run_snakemake_target(
            config_path=config_path,
            cores=int(args.threads),
            target="all",
            snake_args=args.snake_args,
            snakefile=snakefile,
        )
    raise SystemExit(rc)


if __name__ == "__main__":
    main()
