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
RUN_ASR_CHOICES = ("yes", "no")
ALIGNMENT_METHOD_OPTION_MAP: dict[str, tuple[str, ...]] = {
    "1": ("babappalign",),
    "2": ("mafft",),
    "3": ("prank",),
    "4": ("babappalign", "mafft", "prank"),
}
ORTHOGROUP_METHOD_CHOICES = ("rbh", "orthofinder", "rbh_fallback")
ORTHOGROUP_METHOD_LABELS = {
    "rbh": "Reciprocal best hit (RBH) only",
    "orthofinder": "OrthoFinder only",
    "rbh_fallback": "RBH with OrthoFinder fallback/comparison",
}
TRIM_STRATEGY_CHOICES = ("raw", "clipkit", "both")
RECOMBINATION_CHOICES = ("none", "gard", "auto")
GARD_MODE_CHOICES = ("Normal", "Faster")
FORCED_TRIM_STRATEGY = "both"
FORCED_TRIM_STATES = ["raw", "clipkit"]
DEFAULT_TOTAL_THREADS = max(1, os.cpu_count() or 1)
RESUME_STATE_DIRNAME = ".babappasnake"
RESUME_STATE_FILENAME = "resume_state.json"


@dataclass(frozen=True)
class StepSpec:
    rule: str
    description: str
    outputs: tuple[str, ...]


def option_was_provided(raw_args: list[str], option: str) -> bool:
    return any(token == option or token.startswith(f"{option}=") for token in raw_args)


def resume_state_path(outdir: Path) -> Path:
    return outdir / RESUME_STATE_DIRNAME / RESUME_STATE_FILENAME


def load_resume_state(outdir: Path) -> dict[str, object]:
    path = resume_state_path(outdir)
    if not path.exists() or path.stat().st_size == 0:
        return {"completed_steps": {}}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {"completed_steps": {}}
    if not isinstance(payload, dict):
        return {"completed_steps": {}}
    completed = payload.get("completed_steps")
    if not isinstance(completed, dict):
        payload["completed_steps"] = {}
    return payload


def save_resume_state(outdir: Path, payload: dict[str, object]) -> None:
    path = resume_state_path(outdir)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def initialize_resume_state(outdir: Path, config_path: Path, guided_mode: bool) -> None:
    payload = load_resume_state(outdir)
    payload["config_path"] = str(config_path.resolve())
    payload["guided_mode"] = bool(guided_mode)
    payload.setdefault("completed_steps", {})
    save_resume_state(outdir, payload)


def mark_step_completed(outdir: Path, step: StepSpec, status: str) -> None:
    payload = load_resume_state(outdir)
    completed = payload.setdefault("completed_steps", {})
    if not isinstance(completed, dict):
        completed = {}
        payload["completed_steps"] = completed
    completed[step.rule] = {
        "status": status,
        "outputs": list(step.outputs),
    }
    payload["last_completed_step"] = step.rule
    save_resume_state(outdir, payload)


def step_is_recorded_complete(outdir: Path, step: StepSpec) -> bool:
    payload = load_resume_state(outdir)
    completed = payload.get("completed_steps", {})
    return isinstance(completed, dict) and step.rule in completed


def bootstrap_completed_steps_from_outputs(outdir: Path, steps: list[StepSpec]) -> int:
    detected = 0
    for step in steps:
        if step_is_recorded_complete(outdir, step):
            continue
        if step.outputs and all((outdir / rel).exists() for rel in step.outputs):
            mark_step_completed(outdir, step, status="detected_existing_outputs")
            detected += 1
    return detected


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


def parse_trim_strategy(raw: object) -> str:
    strategy = str(raw).strip().lower()
    if strategy not in TRIM_STRATEGY_CHOICES:
        raise SystemExit(
            "Invalid --trim-strategy option. Choose one of: raw, clipkit, both"
        )
    return strategy


def parse_recombination_mode(raw: object) -> str:
    mode = str(raw).strip().lower()
    if mode not in RECOMBINATION_CHOICES:
        raise SystemExit(
            "Invalid --recombination option. Choose one of: none, gard, auto"
        )
    return mode


def effective_recombination_mode(raw: object) -> str:
    mode = parse_recombination_mode(raw)
    return "gard" if mode == "auto" else mode


def derive_trim_strategy(trim_strategy: str | None, use_clipkit: str) -> str:
    if trim_strategy:
        return parse_trim_strategy(trim_strategy)
    return "clipkit" if str(use_clipkit).strip().lower() == "yes" else "raw"


def trim_states_from_strategy(strategy: str) -> list[str]:
    normalized = parse_trim_strategy(strategy)
    if normalized == "raw":
        return ["raw"]
    if normalized == "clipkit":
        return ["clipkit"]
    return ["raw", "clipkit"]


def force_robustness_trim_strategy(requested_strategy: str) -> tuple[str, list[str]]:
    requested = parse_trim_strategy(requested_strategy)
    if requested != FORCED_TRIM_STRATEGY:
        print(
            "[INFO] Trimming strategy is forced to 'both' for robustness "
            "(raw + clipkit summaries will always be generated)."
        )
    return FORCED_TRIM_STRATEGY, list(FORCED_TRIM_STATES)


def enumerate_pathways(methods: list[str], trim_states: list[str]) -> list[str]:
    return [f"{method}_{trim_state}" for method in methods for trim_state in trim_states]


def compute_per_pathway_cores(total_threads: int, pathway_count: int) -> int:
    return max(1, int(total_threads) // max(1, int(pathway_count)))


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


def prompt_orthogroup_method(default: str) -> str:
    print("Orthogroup discovery method (--orthogroup-method)")
    for idx, key in enumerate(ORTHOGROUP_METHOD_CHOICES, start=1):
        mark = " (default)" if key == default else ""
        print(f"  {idx}. {ORTHOGROUP_METHOD_LABELS[key]}{mark}")
    while True:
        value = input("Select option (1/2/3) or press Enter for default: ").strip().lower()
        if not value:
            return default
        if value in {"1", "2", "3"}:
            return ORTHOGROUP_METHOD_CHOICES[int(value) - 1]
        if value in ORTHOGROUP_METHOD_CHOICES:
            return value
        print("Invalid selection.")


def prompt_trim_strategy(default: str) -> str:
    print("Choose trimming strategy")
    options = (
        ("1", "raw", "raw only"),
        ("2", "clipkit", "ClipKIT only"),
        ("3", "both", "both raw and ClipKIT for robustness"),
    )
    for key, strategy, label in options:
        mark = " (default)" if strategy == default else ""
        print(f"  {key}. {label}{mark}")
    while True:
        value = input("Select option (1/2/3) or press Enter for default: ").strip().lower()
        if not value:
            return default
        if value in {"1", "2", "3"}:
            return options[int(value) - 1][1]
        if value in TRIM_STRATEGY_CHOICES:
            return value
        print("Invalid selection.")


def prompt_recombination_mode(default: str) -> str:
    print("Recombination screening (--recombination)")
    options = (
        ("1", "none", "none (default behavior: skip GARD)"),
        ("2", "gard", "run HyPhy GARD screening"),
        ("3", "auto", "auto (currently alias of gard)"),
    )
    for key, mode, label in options:
        mark = " (default)" if mode == default else ""
        print(f"  {key}. {label}{mark}")
    while True:
        value = input("Select option (1/2/3) or press Enter for default: ").strip().lower()
        if not value:
            return default
        if value in {"1", "2", "3"}:
            return options[int(value) - 1][1]
        if value in RECOMBINATION_CHOICES:
            return value
        print("Invalid selection.")


def prompt_gard_mode(default: str) -> str:
    print("HyPhy GARD mode (--gard-mode)")
    for idx, mode in enumerate(GARD_MODE_CHOICES, start=1):
        mark = " (default)" if mode == default else ""
        print(f"  {idx}. {mode}{mark}")
    while True:
        value = input("Select option (1/2) or press Enter for default: ").strip()
        if not value:
            return default
        if value in {"1", "2"}:
            return GARD_MODE_CHOICES[int(value) - 1]
        if value in GARD_MODE_CHOICES:
            return value
        print("Invalid selection.")


def validate_selected_alignment_tools(
    methods: list[str],
    trim_states: list[str],
    resolved_tools: dict[str, str],
) -> None:
    missing: list[str] = []
    if "babappalign" in methods and "babappalign" not in resolved_tools:
        missing.append("babappalign")
    if "mafft" in methods and "mafft" not in resolved_tools:
        missing.append("mafft")
    if "prank" in methods and "prank" not in resolved_tools:
        missing.append("prank")
    if "clipkit" in trim_states and "clipkit" not in resolved_tools:
        missing.append("clipkit")
    if missing:
        raise SystemExit(
            "Selected alignment methods require additional tools not found on PATH: "
            + ", ".join(missing)
        )


def validate_orthogroup_method_tools(method: str, resolved_tools: dict[str, str]) -> None:
    selected = str(method).strip().lower()
    if selected == "rbh":
        missing = [tool for tool in ("blastp", "makeblastdb") if tool not in resolved_tools]
        if missing:
            raise SystemExit(
                "Orthogroup method 'rbh' requires BLAST tools on PATH: "
                + ", ".join(missing)
            )
        return
    if selected == "orthofinder":
        missing = [tool for tool in ("orthofinder", "blastp", "makeblastdb") if tool not in resolved_tools]
        if missing:
            raise SystemExit(
                "Orthogroup method 'orthofinder' requires tools on PATH: "
                + ", ".join(missing)
            )
        return
    if selected == "rbh_fallback":
        missing = [tool for tool in ("blastp", "makeblastdb", "orthofinder") if tool not in resolved_tools]
        if missing:
            raise SystemExit(
                "Orthogroup method 'rbh_fallback' requires RBH tools plus OrthoFinder on PATH: "
                + ", ".join(missing)
            )
        return
    raise SystemExit(f"Unsupported orthogroup method: {method}")


def validate_recombination_tools(mode: str, resolved_tools: dict[str, str]) -> None:
    normalized = effective_recombination_mode(mode)
    if normalized == "gard" and "hyphy" not in resolved_tools:
        raise SystemExit(
            "Recombination mode 'gard' requires HyPhy on PATH, but it was not found."
        )


def maybe_prompt_interactive(args: argparse.Namespace) -> argparse.Namespace:
    if getattr(args, "resume", False):
        return args
    if args.interactive == "no":
        return args
    if not is_tty_interactive():
        return args

    print("\nInteractive babappasnake configuration\n")
    args.prot = prompt_text("Proteomes directory (--prot)", args.prot or "", required=True)
    args.query = prompt_text("Query FASTA (--query)", args.query or "", required=True)
    args.outdir = prompt_text("Output directory (--outdir)", args.outdir, required=True)
    args.orthogroup_method = prompt_orthogroup_method(str(args.orthogroup_method))
    args.alignment_methods = prompt_alignment_method_option(str(args.alignment_methods))
    args.coverage = prompt_float("RBH reciprocal coverage (--coverage)", float(args.coverage))
    args.threads = prompt_int("Total cores (--threads)", int(args.threads))
    default_trim = derive_trim_strategy(args.trim_strategy, args.use_clipkit)
    args.trim_strategy, trim_states = force_robustness_trim_strategy(default_trim)
    args.use_clipkit = "yes"
    print(
        "Trimming strategy: both (raw + ClipKIT branches) [enforced for robustness summaries]."
    )
    args.recombination = prompt_recombination_mode(str(args.recombination))
    if effective_recombination_mode(args.recombination) == "gard":
        args.gard_mode = prompt_gard_mode(str(args.gard_mode))
        args.gard_rate_classes = prompt_int(
            "HyPhy GARD rate classes (--gard-rate-classes)",
            int(args.gard_rate_classes),
        )
    if "clipkit" in trim_states:
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
    args.run_asr = "yes" if prompt_yes_no("Run ASR extraction block after branch-site analysis?", args.run_asr == "yes") else "no"
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
    print(
        f"- orthogroup_method: {args.orthogroup_method} "
        f"({ORTHOGROUP_METHOD_LABELS[str(args.orthogroup_method).strip().lower()]})"
    )
    print(f"- alignment_methods: {format_alignment_option(str(args.alignment_methods))}")
    print(f"- trim_strategy: {args.trim_strategy}")
    print(f"- pathways: {', '.join(enumerate_pathways(parse_alignment_method_option(args.alignment_methods), trim_states))}")
    print(f"- recombination: {args.recombination} (effective: {effective_recombination_mode(args.recombination)})")
    if effective_recombination_mode(args.recombination) == "gard":
        print(f"- gard_mode: {args.gard_mode}")
        print(f"- gard_rate_classes: {args.gard_rate_classes}")
    print(f"- threads: {args.threads}")
    print(f"- coverage: {args.coverage}")
    if args.orthogroup_method == "rbh_fallback":
        print("- orthogroup_selection: compare RBH vs OrthoFinder and keep the higher 1:1 count")
    print(f"- use_clipkit (compat): {args.use_clipkit}")
    if "clipkit" in trim_states:
        print(f"- clipkit_mode_protein: {args.clipkit_mode_protein}")
        print(f"- clipkit_mode_codon: {args.clipkit_mode_codon}")
    print(f"- iqtree_bootstrap: {args.iqtree_bootstrap}")
    print(f"- iqtree_bnni: {args.iqtree_bnni}")
    print(f"- iqtree_model: {args.iqtree_model}")
    print(f"- absrel_branches: {args.absrel_branches}")
    print(f"- meme_branches: {args.meme_branches}")
    print(f"- codeml_codonfreq: {args.codeml_codonfreq}")
    print(f"- run_asr: {args.run_asr}")
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
    p.add_argument(
        "--orthogroup-method",
        choices=list(ORTHOGROUP_METHOD_CHOICES),
        default="rbh",
        help=(
            "Orthogroup discovery backend [default: rbh]. "
            "Choose RBH only, OrthoFinder only, or RBH with OrthoFinder fallback/comparison."
        ),
    )
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
    p.add_argument("--codeml-codonfreq", type=int, default=7, help="codeml CodonFreq value [default: 7]")
    p.add_argument("--run-asr", choices=list(RUN_ASR_CHOICES), default="yes", help="Run the ASR extraction block after branch-site analysis [default: yes]")
    p.add_argument(
        "--recombination",
        choices=list(RECOMBINATION_CHOICES),
        default="none",
        help="Optional recombination screening mode: none, gard, or auto (currently alias of gard) [default: none]",
    )
    p.add_argument(
        "--gard-mode",
        choices=list(GARD_MODE_CHOICES),
        default="Faster",
        help="HyPhy GARD mode when recombination screening is enabled [default: Faster]",
    )
    p.add_argument(
        "--gard-rate-classes",
        type=int,
        default=3,
        help="HyPhy GARD rate classes when recombination screening is enabled [default: 3]",
    )
    p.add_argument(
        "--trim-strategy",
        choices=list(TRIM_STRATEGY_CHOICES),
        default=None,
        help=(
            "Trimming strategy: raw (no ClipKIT), clipkit, or both for robustness. "
            "If omitted, legacy --use-clipkit mapping is used."
        ),
    )
    p.add_argument("--clipkit-mode-protein", default="kpic-smart-gap", help="ClipKIT mode for protein trimming")
    p.add_argument("--clipkit-mode-codon", default="kpic-smart-gap", help="ClipKIT mode for codon trimming")
    p.add_argument("--absrel-p", type=float, default=0.05, help="Leaf-branch significance threshold for aBSREL [default: 0.05]")
    p.add_argument("--absrel-dynamic-start", type=float, default=0.05, help="Dynamic aBSREL start threshold [default: 0.05]")
    p.add_argument("--absrel-dynamic-step", type=float, default=0.01, help="Dynamic aBSREL step size [default: 0.01]")
    p.add_argument("--absrel-dynamic-max", type=float, default=0.2, help="Dynamic aBSREL max threshold [default: 0.2]")
    p.add_argument("--meme-p", type=float, default=0.05, help="MEME site threshold retained in summary [default: 0.05]")
    p.add_argument(
        "--use-clipkit",
        choices=["yes", "no"],
        default="yes",
        help="Legacy compatibility switch mapping to --trim-strategy clipkit/raw when --trim-strategy is not provided [default: yes]",
    )
    p.add_argument("--snake-args", default="", help="Extra raw arguments forwarded to snakemake")
    p.add_argument("--interactive", choices=["yes", "no"], default="yes", help="Prompt for pipeline settings interactively [default: yes]")
    p.add_argument("--guided", choices=["yes", "no"], default="yes", help="Run step-by-step guided execution [default: yes]")
    p.add_argument("--resume", action="store_true", help="Resume an existing run from the saved outdir/config.yaml")
    return p.parse_args()


def load_config_file(config_path: Path) -> dict[str, object]:
    with open(config_path, "r", encoding="utf-8") as fh:
        payload = yaml.safe_load(fh) or {}
    if not isinstance(payload, dict):
        raise SystemExit(f"Invalid config structure in {config_path}")
    return payload


def save_config_file(config_path: Path, payload: dict[str, object]) -> None:
    with open(config_path, "w", encoding="utf-8") as fh:
        yaml.safe_dump(payload, fh, sort_keys=False)


def bool_to_yes_no(value: object) -> str:
    return "yes" if bool(value) else "no"


def alignment_option_from_methods(methods: object) -> str:
    if isinstance(methods, str):
        normalized = tuple(token.strip().lower() for token in methods.split(",") if token.strip())
    elif isinstance(methods, (list, tuple)):
        normalized = tuple(str(token).strip().lower() for token in methods if str(token).strip())
    else:
        normalized = ()
    for option, mapped in ALIGNMENT_METHOD_OPTION_MAP.items():
        if tuple(mapped) == normalized:
            return option
    raise SystemExit(
        "Saved run config contains unsupported alignment methods: "
        + ", ".join(normalized or ("(none)",))
    )


def validate_resume_request(args: argparse.Namespace, raw_args: list[str]) -> Path:
    forbidden = (
        "--prot",
        "--query",
        "--orthogroup-method",
        "--alignment-methods",
        "--coverage",
        "--iqtree-bootstrap",
        "--iqtree-bnni",
        "--iqtree-model",
        "--absrel-branches",
        "--meme-branches",
        "--codeml-codonfreq",
        "--run-asr",
        "--recombination",
        "--gard-mode",
        "--gard-rate-classes",
        "--trim-strategy",
        "--clipkit-mode-protein",
        "--clipkit-mode-codon",
        "--absrel-p",
        "--absrel-dynamic-start",
        "--absrel-dynamic-step",
        "--absrel-dynamic-max",
        "--meme-p",
        "--use-clipkit",
    )
    provided = [opt for opt in forbidden if option_was_provided(raw_args, opt)]
    if provided:
        raise SystemExit(
            "When using --resume, analysis settings are loaded from the saved config. "
            "Do not pass: " + ", ".join(provided)
        )
    if args.cds:
        cds_path = Path(args.cds).expanduser()
        if not cds_path.is_file():
            raise SystemExit(f"--cds must point to an existing FASTA file: {cds_path}")
    outdir = Path(args.outdir)
    config_path = outdir / "config.yaml"
    if not config_path.is_file():
        raise SystemExit(
            f"Cannot resume because no saved config was found at: {config_path}"
        )
    return config_path


def hydrate_args_from_saved_config(
    args: argparse.Namespace,
    config: dict[str, object],
    raw_args: list[str],
) -> argparse.Namespace:
    args.orthogroup_method = str(config.get("orthogroup_method", "rbh")).strip().lower() or "rbh"
    args.alignment_methods = str(
        config.get("alignment_method_option", alignment_option_from_methods(config.get("alignment_methods", [])))
    )
    args.coverage = float(config.get("coverage", args.coverage))
    if not option_was_provided(raw_args, "--threads"):
        args.threads = int(config.get("threads", args.threads))
    args.iqtree_bootstrap = int(config.get("iqtree_bootstrap", args.iqtree_bootstrap))
    args.iqtree_bnni = bool_to_yes_no(config.get("iqtree_bnni", False))
    args.iqtree_model = str(config.get("iqtree_model", args.iqtree_model)).strip() or "MFP"
    args.absrel_branches = str(config.get("absrel_branches", args.absrel_branches)).strip() or "Leaves"
    args.meme_branches = str(config.get("meme_branches", args.meme_branches)).strip() or "Leaves"
    args.codeml_codonfreq = int(config.get("codeml_codonfreq", args.codeml_codonfreq))
    args.run_asr = bool_to_yes_no(config.get("run_asr", True))
    args.recombination = parse_recombination_mode(config.get("recombination", args.recombination))
    args.gard_mode = str(config.get("gard_mode", args.gard_mode)).strip() or "Faster"
    args.gard_rate_classes = int(config.get("gard_rate_classes", args.gard_rate_classes))
    args.trim_strategy = str(config.get("trim_strategy", args.trim_strategy or FORCED_TRIM_STRATEGY)).strip().lower()
    args.use_clipkit = bool_to_yes_no(config.get("use_clipkit", True))
    args.clipkit_mode_protein = str(config.get("clipkit_mode_protein", args.clipkit_mode_protein)).strip() or "kpic-smart-gap"
    args.clipkit_mode_codon = str(config.get("clipkit_mode_codon", args.clipkit_mode_codon)).strip() or "kpic-smart-gap"
    args.absrel_p = float(config.get("absrel_p", args.absrel_p))
    args.absrel_dynamic_start = float(config.get("absrel_dynamic_start", args.absrel_dynamic_start))
    args.absrel_dynamic_step = float(config.get("absrel_dynamic_step", args.absrel_dynamic_step))
    args.absrel_dynamic_max = float(config.get("absrel_dynamic_max", args.absrel_dynamic_max))
    args.meme_p = float(config.get("meme_p", args.meme_p))
    if not option_was_provided(raw_args, "--outgroup"):
        args.outgroup = str(config.get("outgroup_query", args.outgroup)).strip()
    if not option_was_provided(raw_args, "--guided"):
        args.guided = bool_to_yes_no(config.get("guided_mode", args.guided == "yes"))
    if not option_was_provided(raw_args, "--snake-args"):
        args.snake_args = str(config.get("snake_args", args.snake_args or "")).strip()
    args.interactive = "no"
    return args


def stage_resume_overrides(
    args: argparse.Namespace,
    outdir: Path,
    config_path: Path,
    raw_args: list[str],
) -> None:
    config = load_config_file(config_path)
    changed = False
    if args.cds:
        cds_dst_dir = outdir / "user_supplied"
        cds_dst_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(Path(args.cds).expanduser(), cds_dst_dir / "orthogroup_cds.fasta")
    if option_was_provided(raw_args, "--outgroup"):
        config["outgroup_query"] = str(args.outgroup).strip()
        changed = True
    if option_was_provided(raw_args, "--threads"):
        config["threads"] = int(args.threads)
        changed = True
    if option_was_provided(raw_args, "--guided"):
        config["guided_mode"] = args.guided == "yes"
        changed = True
    if option_was_provided(raw_args, "--snake-args"):
        config["snake_args"] = str(args.snake_args).strip()
        changed = True
    methods = config.get("alignment_methods", [])
    trim_states = config.get("trim_states", [])
    try:
        method_count = max(1, len(methods) if isinstance(methods, list) else len(parse_alignment_method_option(str(args.alignment_methods))))
    except Exception:
        method_count = max(1, len(parse_alignment_method_option(str(args.alignment_methods))))
    if isinstance(trim_states, list) and trim_states:
        pathway_count = max(1, method_count * len(trim_states))
    else:
        pathway_count = max(1, method_count * len(trim_states_from_strategy(args.trim_strategy)))
    config["per_method_cores"] = compute_per_method_cores(int(args.threads), method_count)
    config["per_pathway_cores"] = compute_per_pathway_cores(int(args.threads), pathway_count)
    changed = True
    if changed:
        save_config_file(config_path, config)


def prepare_resume_run(
    args: argparse.Namespace,
    raw_args: list[str],
) -> tuple[argparse.Namespace, Path, Path]:
    config_path = validate_resume_request(args, raw_args)
    outdir = Path(args.outdir)
    config = load_config_file(config_path)
    args = hydrate_args_from_saved_config(args, config, raw_args)
    stage_resume_overrides(args, outdir, config_path, raw_args)
    initialize_resume_state(outdir, config_path, args.guided == "yes")
    return args, config_path, outdir


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
    if str(args.orthogroup_method).strip().lower() not in ORTHOGROUP_METHOD_CHOICES:
        raise SystemExit(
            f"Invalid --orthogroup-method option: {args.orthogroup_method}. "
            "Choose one of: rbh, orthofinder, rbh_fallback"
        )
    parse_trim_strategy(derive_trim_strategy(args.trim_strategy, args.use_clipkit))
    parse_recombination_mode(args.recombination)
    if int(args.gard_rate_classes) < 1:
        raise SystemExit("--gard-rate-classes must be >= 1")


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
    trim_strategy: str,
    trim_states: list[str],
    per_method_cores: int,
    per_pathway_cores: int,
) -> Path:
    cfg = {
        "outdir": str(outdir.resolve()),
        "query_fasta": str((outdir / "inputs" / "query.fasta").resolve()),
        "proteomes_dir": str((outdir / "inputs" / "proteomes").resolve()),
        "user_cds": str((outdir / "user_supplied" / "orthogroup_cds.fasta").resolve()),
        "orthogroup_method": str(args.orthogroup_method).strip().lower(),
        "alignment_methods": methods,
        "alignment_method_option": str(args.alignment_methods).strip(),
        "trim_strategy": trim_strategy,
        "trim_states": trim_states,
        "per_method_cores": int(per_method_cores),
        "per_pathway_cores": int(per_pathway_cores),
        "coverage": float(args.coverage),
        "threads": int(args.threads),
        "outgroup_query": str(args.outgroup).strip(),
        "iqtree_bootstrap": int(args.iqtree_bootstrap),
        "iqtree_bnni": args.iqtree_bnni == "yes",
        "iqtree_model": str(args.iqtree_model).strip(),
        "absrel_branches": str(args.absrel_branches).strip() or "Leaves",
        "meme_branches": str(args.meme_branches).strip() or "Leaves",
        "codeml_codonfreq": int(args.codeml_codonfreq),
        "run_asr": args.run_asr == "yes",
        "recombination": effective_recombination_mode(args.recombination),
        "gard_mode": str(args.gard_mode).strip() or "Faster",
        "gard_rate_classes": int(args.gard_rate_classes),
        "use_clipkit": "clipkit" in trim_states,
        "clipkit_mode_protein": str(args.clipkit_mode_protein).strip(),
        "clipkit_mode_codon": str(args.clipkit_mode_codon).strip(),
        "absrel_p": float(args.absrel_p),
        "absrel_dynamic_start": float(args.absrel_dynamic_start),
        "absrel_dynamic_step": float(args.absrel_dynamic_step),
        "absrel_dynamic_max": float(args.absrel_dynamic_max),
        "meme_p": float(args.meme_p),
        "guided_mode": args.guided == "yes",
        "snake_args": str(args.snake_args).strip(),
        "executables": executables,
    }
    config_path = outdir / "config.yaml"
    with open(config_path, "w", encoding="utf-8") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)
    return config_path


def build_snakemake_cmd(
    config_path: Path,
    cores: int,
    target: str,
    snake_args: str,
    snakefile: Path,
    workdir: Path,
) -> list[str]:
    cmd = [
        sys.executable,
        "-m",
        "snakemake",
        "--directory",
        str(workdir.resolve()),
        "--snakefile",
        str(Path(str(snakefile)).resolve()),
        "--configfile",
        str(config_path.resolve()),
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


def unlock_snakemake_run(config_path: Path, snakefile: Path, workdir: Path) -> int:
    cmd = [
        sys.executable,
        "-m",
        "snakemake",
        "--directory",
        str(workdir.resolve()),
        "--snakefile",
        str(Path(str(snakefile)).resolve()),
        "--configfile",
        str(config_path.resolve()),
        "--unlock",
    ]
    result = subprocess.run(cmd, check=False)
    return result.returncode


def run_snakemake_target(
    config_path: Path,
    cores: int,
    target: str,
    snake_args: str,
    snakefile: Path,
    workdir: Path,
) -> int:
    cmd = build_snakemake_cmd(config_path, cores, target, snake_args, snakefile, workdir)
    result = subprocess.run(cmd, check=False)
    return result.returncode


def build_step_plan(
    have_cds: bool,
    methods: list[str],
    trim_states: list[str],
    recombination_mode: str,
    run_asr: bool,
) -> list[StepSpec]:
    pathways = enumerate_pathways(methods, trim_states)

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
        StepSpec(
            "prepare_branch_inputs_all_pathways",
            "Prepare raw/ClipKIT branch inputs for each selected pathway.",
            tuple(
                f"alignments/{method}/{trim_state}/mapped_orthogroup_cds.analysis.fasta"
                for method in methods
                for trim_state in trim_states
            ),
        ),
    ]
    run_gard = effective_recombination_mode(recombination_mode) == "gard"
    if run_gard:
        steps.append(
            StepSpec(
                "gard_all_pathways",
                "Run optional HyPhy GARD recombination screening for each selected pathway.",
                tuple(
                    f"recombination/{method}/{trim_state}/gard/gard_summary.json"
                    for method in methods
                    for trim_state in trim_states
                ),
            )
        )

    later_steps: list[StepSpec] = [
            StepSpec(
                "iqtree_ml_all_pathways",
                "Infer ML phylogeny with IQ-TREE for all selected pathways.",
                tuple(f"tree/{method}/{trim_state}/orthogroup.treefile" for method in methods for trim_state in trim_states),
            ),
            StepSpec(
                "root_iqtree_outgroup_all_pathways",
                "Root pathway-specific trees using the outgroup query.",
                tuple(
                    f"tree/{method}/{trim_state}/orthogroup.rooted.treefile"
                    for method in methods
                    for trim_state in trim_states
                ),
            ),
            StepSpec(
                "hyphy_exploratory_all_pathways",
                "Run HyPhy aBSREL + MEME exploratory tests for all selected pathways.",
                tuple(
                    f"hyphy/{method}/{trim_state}/hyphy_done.json"
                    for method in methods
                    for trim_state in trim_states
                ),
            ),
            StepSpec(
                "parse_foregrounds_all_pathways",
                "Select significant foreground branches for all selected pathways.",
                tuple(
                    f"hyphy/{method}/{trim_state}/significant_foregrounds.tsv"
                    for method in methods
                    for trim_state in trim_states
                ),
            ),
            StepSpec(
                "prepare_foreground_trees_all_pathways",
                "Prepare branch-labeled trees for branch-site codeml (all selected pathways).",
                tuple(
                    f"branchsite/{method}/{trim_state}/foreground_trees.tsv"
                    for method in methods
                    for trim_state in trim_states
                ),
            ),
            StepSpec(
                "branchsite_batch_all_pathways",
                "Run branch-site codeml across selected foregrounds (all pathways).",
                tuple(
                    f"branchsite/{method}/{trim_state}/branchsite_results.tsv"
                    for method in methods
                    for trim_state in trim_states
                ),
            ),
            StepSpec(
                "final_summary_all_pathways",
                "Write pathway-specific episodic selection summaries.",
                tuple(
                    f"summary/{method}/{trim_state}/episodic_selection_summary.txt"
                    for method in methods
                    for trim_state in trim_states
                ),
            ),
            StepSpec(
                "robustness_reports",
                "Write robustness matrix, consensus, narrative, and publication-ready table across pathways.",
                (
                    "summary/robustness_matrix.tsv",
                    "summary/robustness_consensus.tsv",
                    "summary/robustness_narrative.txt",
                    "summary/comparative_reproducibility_summary.txt",
                    "summary/robustness_publication_table.tex",
                ),
            ),
            StepSpec(
                "final_summary_primary_alias",
                "Write legacy top-level final summary alias.",
                ("summary/episodic_selection_summary.txt",),
            ),
            StepSpec(
                "write_run_provenance",
                "Write machine-readable run provenance.",
                ("summary/run_provenance.json",),
            ),
    ]
    if run_asr:
        later_steps[6:6] = [
            StepSpec(
                "codeml_asr_all_pathways",
                "Run codeml ancestral sequence reconstruction for all selected pathways.",
                tuple(
                    f"asr/{method}/{trim_state}/asr_done.json"
                    for method in methods
                    for trim_state in trim_states
                ),
            ),
            StepSpec(
                "extract_selected_branch_ancestors",
                "Map selected branches to parent/child nodes and recover ancestral/descendant sequences with substitutions.",
                (
                    "asr/branch_to_nodes.tsv",
                    "asr/ancestor_sequences_cds.fasta",
                    "asr/ancestor_sequences_aa.fasta",
                    "asr/descendant_sequences_cds.fasta",
                    "asr/descendant_sequences_aa.fasta",
                    "asr/branch_substitutions.tsv",
                    "asr/selected_branch_asr_summary.tsv",
                    "asr/asr_extraction_provenance.json",
                ),
            ),
        ]
        later_steps.insert(
            -1,
            StepSpec(
                "asr_primary_alias",
                "Write legacy top-level ASR completion alias.",
                ("asr/asr_done.json",),
            ),
        )
    steps.extend(later_steps)
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
    if step_is_recorded_complete(outdir, step) and all(existing):
        print("This step is already recorded complete. Auto-skipping.")
        print_step_outputs(outdir, step.outputs)
        return 0
    if all(existing):
        print(
            "All expected outputs already exist, but this step is not yet recorded complete. "
            "Choose 'run' to let Snakemake verify it, or skip manually."
        )
    can_skip = all(existing) or force_can_skip
    if can_skip:
        if force_can_skip:
            print("This step is skippable in the current context.")
    action = prompt_step_action(can_skip=can_skip)
    if action == "stop":
        print("Guided run stopped by user.")
        return 130
    if action == "skip":
        if on_skip is not None:
            on_skip()
        mark_step_completed(outdir, step, status="skipped")
        print("Step skipped (existing outputs retained).")
        print_step_outputs(outdir, step.outputs)
        return 0
    code = run_snakemake_target(
        config_path,
        cores,
        step.rule,
        snake_args,
        snakefile,
        outdir,
    )
    if code != 0:
        return code
    mark_step_completed(outdir, step, status="completed")
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
    print(f"  Included groups with ortholog(s): {len(included)}")
    for species, ortholog in included:
        print(f"    - {species}: {ortholog}")
    print(f"  Omitted groups (no ortholog in selected orthogroup): {len(omitted)}")
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


def maybe_prompt_outgroup_after_cds(
    args: argparse.Namespace,
    config_path: Path,
    have_cds: bool,
) -> None:
    if not have_cds:
        print("No CDS staged yet; deferring outgroup prompt until downstream tree stages are reachable.")
        return

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


def set_outgroup_in_config(config_path: Path, outgroup_query: str) -> None:
    with open(config_path, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh) or {}
    cfg["outgroup_query"] = outgroup_query
    with open(config_path, "w", encoding="utf-8") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)


def stage_unrooted_tree_for_downstream(outdir: Path, methods: list[str], trim_states: list[str]) -> None:
    staged = 0
    for method in methods:
        for trim_state in trim_states:
            src = outdir / "tree" / method / trim_state / "orthogroup.treefile"
            dst = outdir / "tree" / method / trim_state / "orthogroup.rooted.treefile"
            if not src.exists():
                print(f"Cannot stage unrooted tree fallback for {method}_{trim_state}; missing: {src}")
                continue
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
            staged += 1
    print(f"Staged unrooted tree fallback for {staged} pathway(s).")


def run_guided_pipeline(
    args: argparse.Namespace,
    config_path: Path,
    outdir: Path,
    cores: int,
    trim_states: list[str],
    snake_args: str,
    snakefile: Path,
) -> int:
    methods = parse_alignment_method_option(args.alignment_methods)
    pathways = enumerate_pathways(methods, trim_states)
    print("Selected pathways:")
    for pathway in pathways:
        print(f"  - {pathway}")
    rbh_step = StepSpec(
        "rbh_orthogroup",
        f"Build orthogroup proteins using '{args.orthogroup_method}' backend.",
        (
            "orthogroup/orthogroup_proteins.fasta",
            "orthogroup/orthogroup_headers.txt",
            "orthogroup/rbh_summary.tsv",
        ),
    )
    print("Guided mode enabled.")
    if getattr(args, "resume", False):
        detected = bootstrap_completed_steps_from_outputs(outdir, [rbh_step])
        if detected:
            print(f"[INFO] Resume scan detected {detected} previously completed step(s) before guided restart.")
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
        maybe_prompt_outgroup_after_cds(args, config_path, have_cds)

    steps = build_step_plan(have_cds, methods, trim_states, args.recombination, args.run_asr == "yes")
    total_steps = 1 + len(steps)
    if getattr(args, "resume", False):
        detected = bootstrap_completed_steps_from_outputs(outdir, steps)
        if detected:
            print(f"[INFO] Resume scan detected {detected} additional completed step(s).")
    print(f"Planned remaining steps: {len(steps)}")
    for idx, step in enumerate(steps, start=1):
        allow_skip = False
        on_skip = None
        if step.rule == "root_iqtree_outgroup_all_pathways" and not str(args.outgroup).strip():
            allow_skip = True
            on_skip = lambda: stage_unrooted_tree_for_downstream(outdir, methods, trim_states)
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
    raw_args = sys.argv[1:]
    args = parse_args()
    snakefile = ir.files("babappasnake").joinpath("workflow", "Snakefile")

    if args.resume:
        args, config_path, outdir = prepare_resume_run(args, raw_args)
        args.recombination = parse_recombination_mode(args.recombination)
        methods = parse_alignment_method_option(args.alignment_methods)
        effective_trim_strategy, trim_states = force_robustness_trim_strategy(
            derive_trim_strategy(args.trim_strategy, args.use_clipkit)
        )
        args.trim_strategy = effective_trim_strategy
        args.use_clipkit = "yes"
        method_count = max(1, len(methods))
        pathway_count = max(1, len(methods) * len(trim_states))
        per_method_cores = compute_per_method_cores(int(args.threads), method_count)
        per_pathway_cores = compute_per_pathway_cores(int(args.threads), pathway_count)
        print(f"[INFO] Resuming existing run from: {outdir}")
        print(
            f"[INFO] Core distribution: total={args.threads}, selected_methods={method_count}, "
            f"selected_trim_states={len(trim_states)}, pathways={pathway_count}, "
            f"per_method={per_method_cores}, per_pathway={per_pathway_cores}."
        )
        print("[INFO] Selected pathways: " + ", ".join(enumerate_pathways(methods, trim_states)))
        print(
            f"[INFO] Recombination screening mode: {args.recombination} "
            f"(effective: {effective_recombination_mode(args.recombination)})."
        )
        print(f"[INFO] ASR extraction block: {args.run_asr}.")
        print("[INFO] Clearing any stale Snakemake lock for this run directory before resume.")
        unlock_code = unlock_snakemake_run(config_path, snakefile, outdir)
        if unlock_code != 0:
            raise SystemExit(unlock_code)
    else:
        args = maybe_prompt_interactive(args)
        validate_inputs(args)
        args.recombination = parse_recombination_mode(args.recombination)
        recombination_mode = effective_recombination_mode(args.recombination)
        methods = parse_alignment_method_option(args.alignment_methods)
        requested_trim_strategy = derive_trim_strategy(args.trim_strategy, args.use_clipkit)
        effective_trim_strategy, trim_states = force_robustness_trim_strategy(requested_trim_strategy)
        args.trim_strategy = effective_trim_strategy
        args.use_clipkit = "yes"

        method_count = max(1, len(methods))
        pathway_count = max(1, len(methods) * len(trim_states))
        if int(args.threads) < pathway_count:
            required_cores = pathway_count
            print(
                f"[INFO] Increasing --threads from {args.threads} to {required_cores} "
                "so each selected pathway can run in parallel."
            )
            args.threads = required_cores
        per_method_cores = compute_per_method_cores(int(args.threads), method_count)
        per_pathway_cores = compute_per_pathway_cores(int(args.threads), pathway_count)
        print(
            f"[INFO] Core distribution: total={args.threads}, selected_methods={method_count}, "
            f"selected_trim_states={len(trim_states)}, pathways={pathway_count}, "
            f"per_method={per_method_cores}, per_pathway={per_pathway_cores}."
        )
        print("[INFO] Selected pathways: " + ", ".join(enumerate_pathways(methods, trim_states)))
        print(
            f"[INFO] Recombination screening mode: {args.recombination} "
            f"(effective: {recombination_mode})."
        )
        print(f"[INFO] ASR extraction block: {args.run_asr}.")
        if recombination_mode == "gard":
            print(
                f"[INFO] GARD parameters: mode={args.gard_mode}, "
                f"rate_classes={args.gard_rate_classes}."
            )

        resolved, missing = resolve_tools()
        if missing:
            print(format_missing_tools(missing), file=sys.stderr)
            raise SystemExit(2)
        validate_orthogroup_method_tools(args.orthogroup_method, resolved)
        validate_recombination_tools(args.recombination, resolved)
        validate_selected_alignment_tools(methods, trim_states, resolved)

        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        stage_inputs(args, outdir)
        config_path = write_config(
            args=args,
            outdir=outdir,
            executables=resolved,
            methods=methods,
            trim_strategy=effective_trim_strategy,
            trim_states=trim_states,
            per_method_cores=per_method_cores,
            per_pathway_cores=per_pathway_cores,
        )
        initialize_resume_state(outdir, config_path, args.guided == "yes")

    if args.guided == "yes" and is_tty_interactive():
        rc = run_guided_pipeline(
            args=args,
            config_path=config_path,
            outdir=outdir,
            cores=int(args.threads),
            trim_states=trim_states,
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
            workdir=outdir,
        )
    raise SystemExit(rc)


if __name__ == "__main__":
    main()
