from __future__ import annotations

import argparse
import csv
import importlib.resources as ir
import json
import os
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import yaml

from babappasnake.utils import format_missing_tools, resolve_tools


IQTREE_BOOTSTRAP_CHOICES = (1000, 5000, 10000)
CODEML_CODONFREQ_CHOICES = (1, 2, 7)
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
HYPHY_BRANCH_CHOICES = ("Leaves", "Internal", "All")
ALIGNMENT_METHOD_OPTION_MAP: dict[str, tuple[str, ...]] = {
    "1": ("babappalign",),
    "2": ("mafft",),
    "3": ("prank",),
    "4": ("babappalign", "mafft", "prank"),
}
DEFAULT_TOTAL_THREADS = max(1, os.cpu_count() or 1)


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


def prompt_hyphy_branches(label: str, default: str) -> str:
    print(label)
    for idx, option in enumerate(HYPHY_BRANCH_CHOICES, start=1):
        mark = " (default)" if option == default else ""
        print(f"  {idx}. {option}{mark}")
    print("  c. custom text")
    while True:
        value = input("Select option (1/2/3/c) or press Enter for default: ").strip()
        if not value:
            return default
        if value in {"1", "2", "3"}:
            return HYPHY_BRANCH_CHOICES[int(value) - 1]
        if value.lower() == "c":
            custom = input("Enter custom HyPhy branches selector: ").strip()
            if custom:
                return custom
            print("Custom selector cannot be empty.")
            continue
        if value in HYPHY_BRANCH_CHOICES:
            return value
        print("Invalid selection.")


def prompt_codonfreq(default: int) -> int:
    print("codeml CodonFreq setting")
    for idx, option in enumerate(CODEML_CODONFREQ_CHOICES, start=1):
        mark = " (default)" if option == default else ""
        print(f"  {idx}. {option}{mark}")
    print("  c. custom integer")
    while True:
        value = input("Select option (1/2/3/c) or press Enter for default: ").strip().lower()
        if not value:
            return default
        if value in {"1", "2", "3"}:
            return CODEML_CODONFREQ_CHOICES[int(value) - 1]
        if value == "c":
            return prompt_int("Enter custom CodonFreq value", default)
        if value.isdigit():
            return int(value)
        print("Invalid selection.")


def parse_alignment_method_option(raw: object) -> list[str]:
    option = str(raw).strip()
    if option not in ALIGNMENT_METHOD_OPTION_MAP:
        raise SystemExit(
            "Invalid --alignment-methods option. Choose one of: "
            "1 (babappalign), 2 (mafft), 3 (prank), 4 (all three)"
        )
    return list(ALIGNMENT_METHOD_OPTION_MAP[option])


def compute_per_method_cores(total_threads: int, method_count: int) -> int:
    return max(1, int(total_threads) // max(1, int(method_count)))


def format_alignment_option(option: str) -> str:
    methods = ALIGNMENT_METHOD_OPTION_MAP.get(option, ())
    return f"{option} ({', '.join(methods)})" if methods else option


def prompt_alignment_method_option(default: str) -> str:
    print("Alignment method selection (--alignment-methods)")
    labels = {
        "1": "babappalign",
        "2": "mafft",
        "3": "prank",
        "4": "all three",
    }
    for key in ("1", "2", "3", "4"):
        mark = " (default)" if key == default else ""
        print(f"  {key}. {labels[key]}{mark}")
    while True:
        value = input("Select option (1/2/3/4) or press Enter for default: ").strip()
        if not value:
            return default
        if value in ALIGNMENT_METHOD_OPTION_MAP:
            return value
        print("Invalid selection.")


def validate_selected_alignment_tools(methods: list[str], resolved_tools: dict[str, str]) -> None:
    missing: list[str] = []
    if "babappalign" in methods and "babappalign" not in resolved_tools:
        missing.append("babappalign")
    if "mafft" in methods and "mafft" not in resolved_tools:
        missing.append("mafft")
    if "prank" in methods and "prank" not in resolved_tools:
        missing.append("prank")
    if any(method in {"mafft", "prank"} for method in methods) and "pal2nal" not in resolved_tools:
        missing.append("pal2nal")
    if missing:
        raise SystemExit(
            "Selected alignment methods require additional tools not found on PATH: "
            + ", ".join(missing)
        )


def maybe_prompt_interactive(args: argparse.Namespace) -> argparse.Namespace:
    if args.interactive == "no":
        return args
    if not is_tty_interactive():
        return args

    print("\nInteractive babappasnake configuration\n")
    args.prot = prompt_text("Proteomes directory (--prot)", args.prot or "", required=True)
    args.query = prompt_text("Query FASTA (--query)", args.query or "", required=True)
    args.outdir = prompt_text("Output directory (--outdir)", args.outdir, required=True)
    args.alignment_methods = prompt_alignment_method_option(str(args.alignment_methods))
    args.coverage = prompt_float("RBH reciprocal coverage (--coverage)", float(args.coverage))
    args.threads = prompt_int("Total cores (--threads)", int(args.threads))
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
    args.absrel_branches = prompt_hyphy_branches(
        "HyPhy aBSREL branches selector (--absrel-branches)",
        str(args.absrel_branches),
    )
    args.meme_branches = prompt_hyphy_branches(
        "HyPhy MEME branches selector (--meme-branches)",
        str(args.meme_branches),
    )
    args.codeml_codonfreq = prompt_codonfreq(int(args.codeml_codonfreq))
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
    print(f"- outdir: {args.outdir}")
    print(f"- alignment_methods: {format_alignment_option(str(args.alignment_methods))}")
    print(f"- threads: {args.threads}")
    print(f"- coverage: {args.coverage}")
    print(f"- use_clipkit: {args.use_clipkit}")
    if args.use_clipkit == "yes":
        print(f"- clipkit_mode_protein: {args.clipkit_mode_protein}")
        print(f"- clipkit_mode_codon: {args.clipkit_mode_codon}")
    print(f"- iqtree_bootstrap: {args.iqtree_bootstrap}")
    print(f"- iqtree_bnni: {args.iqtree_bnni}")
    print(f"- iqtree_model: {args.iqtree_model}")
    print(f"- absrel_branches: {args.absrel_branches}")
    print(f"- meme_branches: {args.meme_branches}")
    print(f"- codeml_codonfreq: {args.codeml_codonfreq}")
    print("- cds: prompted after orthogroup definition")
    print("- outgroup: prompted after CDS step (optional)")
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
    p.add_argument(
        "--alignment-methods",
        default="4",
        choices=["1", "2", "3", "4"],
        help="Alignment engine selection: 1=babappalign, 2=mafft, 3=prank, 4=all three [default: 4]",
    )
    p.add_argument("--coverage", type=float, default=0.70, help="RBH minimum reciprocal coverage [default: 0.70]")
    p.add_argument(
        "--threads",
        type=int,
        default=DEFAULT_TOTAL_THREADS,
        help=f"Total Snakemake cores (equally split across selected MSA pathways) [default: {DEFAULT_TOTAL_THREADS}]",
    )
    p.add_argument("--outgroup", default="", help="Outgroup label query used to root the IQ-TREE output (case-insensitive substring match)")
    p.add_argument("--iqtree-bootstrap", type=int, default=1000, help="IQ-TREE ultrafast bootstrap replicates [default: 1000]")
    p.add_argument("--iqtree-bnni", choices=["yes", "no"], default="no", help="Enable IQ-TREE -bnni [default: no]")
    p.add_argument("--iqtree-model", default="MFP", help="IQ-TREE model string [default: MFP]")
    p.add_argument("--absrel-branches", default="Leaves", help="HyPhy aBSREL branches selector [default: Leaves]")
    p.add_argument("--meme-branches", default="Leaves", help="HyPhy MEME branches selector [default: Leaves]")
    p.add_argument("--codeml-codonfreq", type=int, default=2, help="codeml CodonFreq value [default: 2]")
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
    parse_alignment_method_option(args.alignment_methods)


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


def write_config(
    args: argparse.Namespace,
    outdir: Path,
    executables: dict[str, str],
    methods: list[str],
    per_method_cores: int,
) -> Path:
    cfg = {
        "outdir": str(outdir.resolve()),
        "query_fasta": str((outdir / "inputs" / "query.fasta").resolve()),
        "proteomes_dir": str((outdir / "inputs" / "proteomes").resolve()),
        "user_cds": str((outdir / "user_supplied" / "orthogroup_cds.fasta").resolve()),
        "alignment_methods": methods,
        "per_method_cores": int(per_method_cores),
        "coverage": float(args.coverage),
        "threads": int(args.threads),
        "outgroup_query": str(args.outgroup).strip(),
        "iqtree_bootstrap": int(args.iqtree_bootstrap),
        "iqtree_bnni": args.iqtree_bnni == "yes",
        "iqtree_model": str(args.iqtree_model).strip(),
        "absrel_branches": str(args.absrel_branches).strip() or "Leaves",
        "meme_branches": str(args.meme_branches).strip() or "Leaves",
        "codeml_codonfreq": int(args.codeml_codonfreq),
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


def build_step_plan(have_cds: bool, use_clipkit: bool, methods: list[str]) -> list[StepSpec]:
    if not have_cds:
        return [
            StepSpec(
                "waiting_note",
                "Write waiting note for later CDS submission.",
                ("orthogroup/WAITING_FOR_CDS.txt",),
            ),
        ]

    steps: list[StepSpec] = [
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
            "align_proteins_all_methods",
            "Generate protein alignments for selected methods.",
            tuple(f"alignments/{method}/orthogroup_proteins.protein.aln.fasta" for method in methods),
        ),
        StepSpec(
            "align_cds_all_methods",
            "Generate codon alignments for selected methods.",
            tuple(f"alignments/{method}/mapped_orthogroup_cds.codon.aln.fasta" for method in methods),
        ),
    ]
    if use_clipkit:
        steps.extend(
            [
                StepSpec(
                    "trim_protein_alignment_all_methods",
                    "Trim protein alignments with ClipKIT.",
                    tuple(f"trimmed/{method}/orthogroup_proteins.clipkit.fasta" for method in methods),
                ),
                StepSpec(
                    "trim_codon_alignment_all_methods",
                    "Trim codon alignments with ClipKIT.",
                    tuple(f"trimmed/{method}/mapped_orthogroup_cds.clipkit.fasta" for method in methods),
                ),
                StepSpec(
                    "strip_terminal_stop_codon_all_methods",
                    "Remove terminal stop codons from trimmed codon alignments.",
                    tuple(f"trimmed/{method}/mapped_orthogroup_cds.clipkit.nostop.fasta" for method in methods),
                ),
            ]
        )
    steps.extend(
        [
            StepSpec(
                "iqtree_ml_all_methods",
                "Infer ML phylogeny with IQ-TREE for selected methods.",
                tuple(f"tree/{method}/orthogroup.treefile" for method in methods),
            ),
            StepSpec(
                "root_iqtree_outgroup_all_methods",
                "Root method-specific trees using the outgroup query.",
                tuple(f"tree/{method}/orthogroup.rooted.treefile" for method in methods),
            ),
            StepSpec(
                "hyphy_exploratory_all_methods",
                "Run HyPhy aBSREL + MEME exploratory tests for selected methods.",
                tuple(f"hyphy/{method}/hyphy_done.json" for method in methods),
            ),
            StepSpec(
                "parse_foregrounds_all_methods",
                "Select significant foreground branches for selected methods.",
                tuple(f"hyphy/{method}/significant_foregrounds.tsv" for method in methods),
            ),
            StepSpec(
                "prepare_foreground_trees_all_methods",
                "Prepare branch-labeled trees for branch-site codeml (all methods).",
                tuple(f"branchsite/{method}/foreground_trees.tsv" for method in methods),
            ),
            StepSpec(
                "branchsite_batch_all_methods",
                "Run branch-site codeml across selected foregrounds (all methods).",
                tuple(f"branchsite/{method}/branchsite_results.tsv" for method in methods),
            ),
            StepSpec(
                "codeml_asr_all_methods",
                "Run codeml ancestral sequence reconstruction for selected methods.",
                tuple(f"asr/{method}/asr_done.json" for method in methods),
            ),
            StepSpec(
                "final_summary_all_methods",
                "Write method-specific episodic selection summaries.",
                tuple(f"summary/{method}/episodic_selection_summary.txt" for method in methods),
            ),
            StepSpec(
                "final_summary_primary_alias",
                "Write legacy top-level final summary alias.",
                ("summary/episodic_selection_summary.txt",),
            ),
            StepSpec(
                "asr_primary_alias",
                "Write legacy top-level ASR completion alias.",
                ("asr/asr_done.json",),
            ),
            StepSpec(
                "compare_alignment_methods",
                "Write comparative reproducibility summary across alignment methods.",
                ("summary/comparative_reproducibility_summary.txt",),
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


def run_guided_step(
    config_path: Path,
    outdir: Path,
    cores: int,
    snake_args: str,
    snakefile: Path,
    step: StepSpec,
    index: int,
    total: int | None,
    force_can_skip: bool = False,
    on_skip: object | None = None,
) -> int:
    if total is None:
        print(f"\nStep {index}: {step.rule}")
    else:
        print(f"\nStep {index}/{total}: {step.rule}")
    print(step.description)
    existing = tuple((outdir / rel).exists() for rel in step.outputs)
    can_skip = all(existing) or force_can_skip
    if can_skip:
        if all(existing):
            print("All expected outputs for this step already exist.")
        elif force_can_skip:
            print("This step is skippable in the current context.")
    action = prompt_step_action(can_skip=can_skip)
    if action == "stop":
        print("Guided run stopped by user.")
        return 130
    if action == "skip":
        if on_skip is not None:
            on_skip()
        print("Step skipped (existing outputs retained).")
        print_step_outputs(outdir, step.outputs)
        return 0
    code = run_snakemake_target(config_path, cores, step.rule, snake_args, snakefile)
    if code != 0:
        return code
    print_step_outputs(outdir, step.outputs)
    return 0


def print_orthogroup_membership(outdir: Path) -> None:
    summary = outdir / "orthogroup" / "rbh_summary.tsv"
    if not summary.exists():
        print("Orthogroup membership report unavailable (rbh_summary.tsv missing).")
        return

    included: list[tuple[str, str]] = []
    omitted: list[str] = []
    try:
        with open(summary, "r", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                species = str(row.get("species", "")).strip() or "(unknown)"
                ortholog = str(row.get("ortholog", "")).strip()
                if ortholog and ortholog != "NA":
                    included.append((species, ortholog))
                else:
                    omitted.append(species)
    except Exception as exc:
        print(f"Orthogroup membership report unavailable (parse error: {exc}).")
        return

    print("\nOrthogroup membership summary:")
    print(f"  Included groups with RBH orthologs: {len(included)}")
    for species, ortholog in included:
        print(f"    - {species}: {ortholog}")
    print(f"  Omitted groups (no RBH hit at current thresholds): {len(omitted)}")
    for species in omitted:
        print(f"    - {species}")


def prompt_for_cds_after_rbh(args: argparse.Namespace, outdir: Path) -> bool:
    user_cds = outdir / "user_supplied" / "orthogroup_cds.fasta"
    existing = user_cds.exists() and user_cds.stat().st_size > 0
    default = str(args.cds) if args.cds else (str(user_cds) if existing else "")
    print("")
    print("CDS input step (required for codon/tree/hyphy/codeml downstream).")
    while True:
        entered = prompt_text(
            "Provide CDS FASTA path now (or press Enter to continue without CDS)",
            default,
            required=False,
        ).strip()
        if not entered:
            return user_cds.exists() and user_cds.stat().st_size > 0

        candidate = Path(entered).expanduser()
        if not candidate.is_file():
            print(f"CDS file not found: {candidate}")
            default = ""
            continue

        user_cds.parent.mkdir(parents=True, exist_ok=True)
        if candidate.resolve() != user_cds.resolve():
            shutil.copy2(candidate, user_cds)
        print(f"Staged CDS file: {user_cds}")
        return True


def set_outgroup_in_config(config_path: Path, outgroup_query: str) -> None:
    with open(config_path, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh) or {}
    cfg["outgroup_query"] = outgroup_query
    with open(config_path, "w", encoding="utf-8") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)


def stage_unrooted_tree_for_downstream(outdir: Path, methods: list[str]) -> None:
    staged = 0
    for method in methods:
        src = outdir / "tree" / method / "orthogroup.treefile"
        dst = outdir / "tree" / method / "orthogroup.rooted.treefile"
        if not src.exists():
            print(f"Cannot stage unrooted tree fallback for {method}; missing: {src}")
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        staged += 1
    print(f"Staged unrooted tree fallback for {staged} method(s).")


def run_guided_pipeline(
    args: argparse.Namespace,
    config_path: Path,
    outdir: Path,
    cores: int,
    use_clipkit: bool,
    snake_args: str,
    snakefile: Path,
) -> int:
    methods = parse_alignment_method_option(args.alignment_methods)
    rbh_step = StepSpec(
        "rbh_orthogroup",
        "Build orthogroup proteins and RBH mapping.",
        (
            "orthogroup/orthogroup_proteins.fasta",
            "orthogroup/orthogroup_headers.txt",
            "orthogroup/rbh_summary.tsv",
        ),
    )
    print("Guided mode enabled.")
    code = run_guided_step(
        config_path=config_path,
        outdir=outdir,
        cores=cores,
        snake_args=snake_args,
        snakefile=snakefile,
        step=rbh_step,
        index=1,
        total=None,
    )
    if code != 0:
        return code

    print_orthogroup_membership(outdir)

    staged_cds = outdir / "user_supplied" / "orthogroup_cds.fasta"
    have_cds = staged_cds.exists() and staged_cds.stat().st_size > 0
    if is_tty_interactive():
        have_cds = prompt_for_cds_after_rbh(args, outdir)
        outgroup_prompt = prompt_text(
            "Optional outgroup query for rooting (press Enter to skip)",
            str(args.outgroup or ""),
            required=False,
        ).strip()
        args.outgroup = outgroup_prompt
        set_outgroup_in_config(config_path, args.outgroup)
        if args.outgroup:
            print(f"Using outgroup query: {args.outgroup}")
        else:
            print("No outgroup query provided; tree will be left unrooted for downstream use.")

    steps = build_step_plan(have_cds, use_clipkit, methods)
    total_steps = 1 + len(steps)
    print(f"Planned remaining steps: {len(steps)}")
    for idx, step in enumerate(steps, start=1):
        allow_skip = False
        on_skip = None
        if step.rule == "root_iqtree_outgroup_all_methods" and not str(args.outgroup).strip():
            allow_skip = True
            on_skip = lambda: stage_unrooted_tree_for_downstream(outdir, methods)
            print("No outgroup provided: running this step or skipping it will both use the unrooted tree downstream.")
        code = run_guided_step(
            config_path=config_path,
            outdir=outdir,
            cores=cores,
            snake_args=snake_args,
            snakefile=snakefile,
            step=step,
            index=idx + 1,
            total=total_steps,
            force_can_skip=allow_skip,
            on_skip=on_skip,
        )
        if code != 0:
            return code
    return 0


def main() -> None:
    args = parse_args()
    args = maybe_prompt_interactive(args)
    validate_inputs(args)
    methods = parse_alignment_method_option(args.alignment_methods)
    method_count = max(1, len(methods))
    if int(args.threads) < method_count:
        required_cores = method_count
        print(
            f"[INFO] Increasing --threads from {args.threads} to {required_cores} "
            "so each selected MSA pathway can run in parallel."
        )
        args.threads = required_cores
    per_method_cores = compute_per_method_cores(int(args.threads), method_count)
    print(
        f"[INFO] Core distribution: total={args.threads}, selected_methods={method_count}, "
        f"per_method={per_method_cores}."
    )

    resolved, missing = resolve_tools()
    if missing:
        print(format_missing_tools(missing), file=sys.stderr)
        raise SystemExit(2)
    validate_selected_alignment_tools(methods, resolved)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stage_inputs(args, outdir)
    config_path = write_config(
        args=args,
        outdir=outdir,
        executables=resolved,
        methods=methods,
        per_method_cores=per_method_cores,
    )
    snakefile = ir.files("babappasnake").joinpath("workflow", "Snakefile")

    if args.guided == "yes" and is_tty_interactive():
        rc = run_guided_pipeline(
            args=args,
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
