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
ORTHOGROUP_METHOD_CHOICES = ("orthofinder",)
ORTHOLOGY_MODE_CHOICES = ("strict", "representative", "paralog")
ORTHOGROUP_SOURCE_CHOICES = ("orthofinder", "external")
TRIM_STRATEGY_CHOICES = ("raw", "clipkit", "both")
RECOMBINATION_CHOICES = ("none", "gard", "auto")
GARD_MODE_CHOICES = ("Normal", "Faster")
TREE_MODE_CHOICES = ("iqtree", "user")
FORCED_TRIM_STRATEGY = "both"
FORCED_TRIM_STATES = ["raw", "clipkit"]
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


def parse_tree_mode(raw: object) -> str:
    mode = str(raw or "iqtree").strip().lower()
    if mode not in TREE_MODE_CHOICES:
        raise SystemExit("Invalid --tree-mode option. Choose one of: iqtree, user")
    return mode


def tree_mode_from_args(args: argparse.Namespace) -> str:
    if str(getattr(args, "user_tree", "") or "").strip():
        return "user"
    return parse_tree_mode(getattr(args, "tree_mode", "iqtree"))


def parse_yes_no_bool(raw: object, default: bool = True) -> bool:
    text = str(raw).strip().lower()
    if text in {"yes", "y", "true", "1", "on"}:
        return True
    if text in {"no", "n", "false", "0", "off"}:
        return False
    return default


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


def alignment_option_from_methods(methods: list[str]) -> str:
    normalized = tuple(str(method).strip().lower() for method in methods if str(method).strip())
    for option, option_methods in ALIGNMENT_METHOD_OPTION_MAP.items():
        if normalized == option_methods:
            return option
    return "4"


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
    print("Orthogroup discovery method: OrthoFinder")
    return "orthofinder"


def prompt_orthology_mode(default: str) -> str:
    print("OrthoFinder orthology/paralogy mode (--orthology-mode)")
    labels = {
        "strict": "strict single-copy orthologs only",
        "representative": "best representative from multi-copy species",
        "paralog": "retain all paralog copies in the selected orthogroup",
    }
    for idx, key in enumerate(ORTHOLOGY_MODE_CHOICES, start=1):
        mark = " (default)" if key == default else ""
        print(f"  {idx}. {labels[key]}{mark}")
    while True:
        value = input("Select option (1/2/3) or press Enter for default: ").strip().lower()
        if not value:
            return default
        if value in {"1", "2", "3"}:
            return ORTHOLOGY_MODE_CHOICES[int(value) - 1]
        if value in ORTHOLOGY_MODE_CHOICES:
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
    if selected == "orthofinder":
        missing = [tool for tool in ("orthofinder", "blastp", "makeblastdb") if tool not in resolved_tools]
        if missing:
            raise SystemExit(
                "Orthogroup method 'orthofinder' requires tools on PATH: "
                + ", ".join(missing)
        )
        return
    raise SystemExit("Unsupported orthogroup method. Supported: orthofinder")


def orthogroup_source_from_args(args: argparse.Namespace) -> str:
    return "external" if str(getattr(args, "orthogroup_proteins", "") or "").strip() else "orthofinder"


def validate_recombination_tools(mode: str, resolved_tools: dict[str, str]) -> None:
    normalized = effective_recombination_mode(mode)
    if normalized == "gard" and "hyphy" not in resolved_tools:
        raise SystemExit(
            "Recombination mode 'gard' requires HyPhy on PATH, but it was not found."
        )


def filter_missing_tools_for_run(
    missing: list[tuple[object, tuple[str, ...]]],
    tree_mode: str,
) -> list[tuple[object, tuple[str, ...]]]:
    skipped = {"iqtree"} if parse_tree_mode(tree_mode) == "user" else set()
    return [(spec, choices) for spec, choices in missing if getattr(spec, "key", "") not in skipped]


def maybe_prompt_interactive(args: argparse.Namespace) -> argparse.Namespace:
    if args.interactive == "no":
        return args
    if not is_tty_interactive():
        return args

    print("\nInteractive babappasnake configuration\n")
    args.orthogroup_proteins = prompt_text(
        "Externally curated orthogroup protein FASTA (--orthogroup-proteins, optional)",
        args.orthogroup_proteins or "",
        required=False,
    )
    external_orthogroup = bool(str(args.orthogroup_proteins or "").strip())
    args.prot = prompt_text("Proteomes directory (--prot)", args.prot or "", required=not external_orthogroup)
    args.query = prompt_text("Query FASTA (--query)", args.query or "", required=not external_orthogroup)
    args.outdir = prompt_text("Output directory (--outdir)", args.outdir, required=True)
    args.orthogroup_method = "orthofinder"
    if external_orthogroup:
        print("Orthogroup source: external curated protein FASTA; OrthoFinder will be skipped.")
    else:
        print("Orthogroup backend is fixed in interactive mode: OrthoFinder.")
        args.orthology_mode = prompt_orthology_mode(str(args.orthology_mode))
    args.alignment_methods = prompt_alignment_method_option(str(args.alignment_methods))
    args.coverage = prompt_float("Minimum query coverage for OrthoFinder mapping (--coverage)", float(args.coverage))
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
    use_user_tree = prompt_yes_no(
        "Use a user-supplied tree instead of building one with IQ-TREE?",
        tree_mode_from_args(args) == "user",
    )
    if use_user_tree:
        args.tree_mode = "user"
        args.user_tree = prompt_text(
            "User-supplied Newick tree path (--tree)",
            args.user_tree or "",
            required=True,
        )
        print("Tree source: user-supplied Newick tree; IQ-TREE tree inference will be skipped.")
    else:
        args.tree_mode = "iqtree"
        args.user_tree = ""
        args.iqtree_bootstrap = prompt_bootstrap(int(args.iqtree_bootstrap))
        args.iqtree_bnni = "yes" if prompt_yes_no("Use IQ-TREE -bnni?", args.iqtree_bnni == "yes") else "no"
        args.iqtree_model = prompt_text("IQ-TREE model string (--iqtree-model)", args.iqtree_model, required=True)
    args.tree_choice_confirmed = True
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
    print(f"- orthogroup_source: {orthogroup_source_from_args(args)}")
    print(f"- orthogroup_method: {args.orthogroup_method}")
    print(f"- orthology_mode: {args.orthology_mode}")
    print(f"- alignment_methods: {format_alignment_option(str(args.alignment_methods))}")
    print(f"- trim_strategy: {args.trim_strategy}")
    print(f"- pathways: {', '.join(enumerate_pathways(parse_alignment_method_option(args.alignment_methods), trim_states))}")
    print(f"- recombination: {args.recombination} (effective: {effective_recombination_mode(args.recombination)})")
    if effective_recombination_mode(args.recombination) == "gard":
        print(f"- gard_mode: {args.gard_mode}")
        print(f"- gard_rate_classes: {args.gard_rate_classes}")
    print(f"- threads: {args.threads}")
    print(f"- coverage: {args.coverage}")
    print(f"- use_clipkit (compat): {args.use_clipkit}")
    if "clipkit" in trim_states:
        print(f"- clipkit_mode_protein: {args.clipkit_mode_protein}")
        print(f"- clipkit_mode_codon: {args.clipkit_mode_codon}")
    print(f"- tree_mode: {args.tree_mode}")
    if args.tree_mode == "user":
        print(f"- tree: {args.user_tree}")
    else:
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
    p.add_argument(
        "--orthogroup-proteins",
        default=None,
        help=(
            "Externally curated orthogroup protein FASTA. When provided, "
            "BABAPPASNAKE skips OrthoFinder and starts downstream analysis from this FASTA."
        ),
    )
    p.add_argument(
        "--orthogroup-method",
        choices=list(ORTHOGROUP_METHOD_CHOICES),
        default="orthofinder",
        help="Orthogroup discovery backend [default: orthofinder].",
    )
    p.add_argument(
        "--orthology-mode",
        choices=list(ORTHOLOGY_MODE_CHOICES),
        default="representative",
        help=(
            "OrthoFinder extraction mode: strict=single-copy only, "
            "representative=best copy from multi-copy species, paralog=retain all copies "
            "[default: representative]"
        ),
    )
    p.add_argument("--outdir", default="babappasnake_run", help="Output directory")
    p.add_argument(
        "--alignment-methods",
        default="4",
        choices=["1", "2", "3", "4"],
        help="Alignment engine selection: 1=babappalign, 2=mafft, 3=prank, 4=all three [default: 4]",
    )
    p.add_argument("--coverage", type=float, default=0.70, help="Minimum query coverage for OrthoFinder mapping [default: 0.70]")
    p.add_argument(
        "--threads",
        type=int,
        default=DEFAULT_TOTAL_THREADS,
        help=f"Total Snakemake cores (equally split across selected MSA pathways) [default: {DEFAULT_TOTAL_THREADS}]",
    )
    p.add_argument("--outgroup", default="", help="Outgroup label query used to root the IQ-TREE output (case-insensitive substring match)")
    p.add_argument(
        "--tree-mode",
        choices=list(TREE_MODE_CHOICES),
        default="iqtree",
        help="Tree source: iqtree builds pathway trees, user reuses a supplied Newick tree for all pathways [default: iqtree]",
    )
    p.add_argument(
        "--tree",
        "--user-tree",
        dest="user_tree",
        default=None,
        help="User-supplied Newick tree used for all pathways when --tree-mode user is selected",
    )
    p.add_argument("--iqtree-bootstrap", type=int, default=1000, help="IQ-TREE ultrafast bootstrap replicates [default: 1000]")
    p.add_argument("--iqtree-bnni", choices=["yes", "no"], default="no", help="Enable IQ-TREE -bnni [default: no]")
    p.add_argument("--iqtree-model", default="MFP", help="IQ-TREE model string [default: MFP]")
    p.add_argument("--absrel-branches", default="Leaves", help="HyPhy aBSREL branches selector [default: Leaves]")
    p.add_argument("--meme-branches", default="Leaves", help="HyPhy MEME branches selector [default: Leaves]")
    p.add_argument("--codeml-codonfreq", type=int, default=7, help="codeml CodonFreq value [default: 7]")
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
    p.add_argument("--run-asr", choices=["yes", "no"], default="yes", help="Run codeml ASR and selected-branch extraction [default: yes]")
    p.add_argument(
        "--use-clipkit",
        choices=["yes", "no"],
        default="yes",
        help="Legacy compatibility switch mapping to --trim-strategy clipkit/raw when --trim-strategy is not provided [default: yes]",
    )
    p.add_argument("--snake-args", default="", help="Extra raw arguments forwarded to snakemake")
    p.add_argument("--interactive", choices=["yes", "no"], default="yes", help="Prompt for pipeline settings interactively [default: yes]")
    p.add_argument("--guided", choices=["yes", "no"], default="yes", help="Run step-by-step guided execution [default: yes]")
    p.add_argument("--resume", action="store_true", help="Resume an existing run directory")
    return p.parse_args()


def validate_inputs(args: argparse.Namespace) -> None:
    missing = []
    external_orthogroup = orthogroup_source_from_args(args) == "external"
    if not external_orthogroup and not args.prot:
        missing.append("--prot")
    if not external_orthogroup and not args.query:
        missing.append("--query")
    if missing:
        raise SystemExit(f"Missing required arguments: {', '.join(missing)}")

    if args.orthogroup_proteins:
        orthogroup_path = Path(args.orthogroup_proteins)
        if not orthogroup_path.is_file():
            raise SystemExit(f"--orthogroup-proteins must point to an existing FASTA file: {orthogroup_path}")
    if args.prot:
        prot_path = Path(args.prot)
        if not prot_path.is_dir():
            raise SystemExit(f"--prot must be an existing directory: {prot_path}")
    if args.query:
        query_path = Path(args.query)
        if not query_path.is_file():
            raise SystemExit(f"--query must be an existing FASTA file: {query_path}")
    if args.cds:
        cds_path = Path(args.cds)
        if not cds_path.is_file():
            raise SystemExit(f"--cds must point to an existing FASTA file: {cds_path}")
    tree_mode = tree_mode_from_args(args)
    args.tree_mode = tree_mode
    if tree_mode == "user":
        user_tree = Path(str(getattr(args, "user_tree", "") or "")).expanduser()
        if not str(getattr(args, "user_tree", "") or "").strip():
            raise SystemExit("--tree-mode user requires --tree / --user-tree pointing to an existing Newick tree.")
        if not user_tree.is_file():
            raise SystemExit(f"--tree must point to an existing Newick tree file: {user_tree}")
    parse_alignment_method_option(args.alignment_methods)
    if str(args.orthogroup_method).strip().lower() not in ORTHOGROUP_METHOD_CHOICES:
        raise SystemExit(
            f"Invalid --orthogroup-method option: {args.orthogroup_method}. "
            "Choose: orthofinder"
        )
    if str(args.orthology_mode).strip().lower() not in ORTHOLOGY_MODE_CHOICES:
        raise SystemExit(
            f"Invalid --orthology-mode option: {args.orthology_mode}. "
            "Choose one of: strict, representative, paralog"
        )
    parse_trim_strategy(derive_trim_strategy(args.trim_strategy, args.use_clipkit))
    parse_recombination_mode(args.recombination)
    if int(args.gard_rate_classes) < 1:
        raise SystemExit("--gard-rate-classes must be >= 1")


def stage_inputs(args: argparse.Namespace, outdir: Path) -> None:
    staged = outdir / "inputs"
    staged.mkdir(parents=True, exist_ok=True)
    if args.query:
        shutil.copy2(args.query, staged / "query.fasta")
    if args.prot:
        prot_dst = staged / "proteomes"
        prot_dst.mkdir(exist_ok=True)
        for item in sorted(Path(args.prot).iterdir()):
            if item.is_file():
                shutil.copy2(item, prot_dst / item.name)
    if args.orthogroup_proteins:
        external_dst_dir = outdir / "user_supplied"
        external_dst_dir.mkdir(exist_ok=True)
        shutil.copy2(args.orthogroup_proteins, external_dst_dir / "external_orthogroup_proteins.fasta")
    if args.cds:
        cds_dst_dir = outdir / "user_supplied"
        cds_dst_dir.mkdir(exist_ok=True)
        shutil.copy2(args.cds, cds_dst_dir / "orthogroup_cds.fasta")
    if getattr(args, "user_tree", None):
        tree_dst_dir = outdir / "user_supplied"
        tree_dst_dir.mkdir(exist_ok=True)
        shutil.copy2(Path(str(args.user_tree)).expanduser(), tree_dst_dir / "user_tree.nwk")


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
    query_path = outdir / "inputs" / "query.fasta"
    proteomes_path = outdir / "inputs" / "proteomes"
    external_orthogroup_path = outdir / "user_supplied" / "external_orthogroup_proteins.fasta"
    user_tree_path = outdir / "user_supplied" / "user_tree.nwk"
    orthogroup_source = "external" if external_orthogroup_path.exists() else "orthofinder"
    tree_mode = "user" if user_tree_path.exists() else parse_tree_mode(getattr(args, "tree_mode", "iqtree"))
    cfg = {
        "outdir": str(outdir.resolve()),
        "query_fasta": str(query_path.resolve()) if query_path.exists() else "",
        "proteomes_dir": str(proteomes_path.resolve()) if proteomes_path.exists() else "",
        "user_cds": str((outdir / "user_supplied" / "orthogroup_cds.fasta").resolve()),
        "orthogroup_source": orthogroup_source,
        "external_orthogroup_proteins": (
            str(external_orthogroup_path.resolve()) if external_orthogroup_path.exists() else ""
        ),
        "orthogroup_method": "orthofinder",
        "orthology_mode": str(getattr(args, "orthology_mode", "representative")).strip().lower(),
        "alignment_methods": methods,
        "alignment_method_option": str(args.alignment_methods),
        "trim_strategy": trim_strategy,
        "trim_states": trim_states,
        "per_method_cores": int(per_method_cores),
        "per_pathway_cores": int(per_pathway_cores),
        "coverage": float(args.coverage),
        "threads": int(args.threads),
        "outgroup_query": str(args.outgroup).strip(),
        "tree_mode": tree_mode,
        "user_tree": str(user_tree_path.resolve()) if user_tree_path.exists() else "",
        "iqtree_bootstrap": int(args.iqtree_bootstrap),
        "iqtree_bnni": args.iqtree_bnni == "yes",
        "iqtree_model": str(args.iqtree_model).strip(),
        "absrel_branches": str(args.absrel_branches).strip() or "Leaves",
        "meme_branches": str(args.meme_branches).strip() or "Leaves",
        "codeml_codonfreq": int(args.codeml_codonfreq),
        "run_asr": parse_yes_no_bool(getattr(args, "run_asr", "yes"), True),
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
        "guided_mode": str(getattr(args, "guided", "yes")).strip().lower() == "yes",
        "snake_args": str(getattr(args, "snake_args", "") or ""),
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
    workdir: Path | None = None,
) -> list[str]:
    run_dir = Path(workdir) if workdir is not None else config_path.parent
    cmd = [
        sys.executable,
        "-m",
        "snakemake",
        "--snakefile",
        str(snakefile),
        "--directory",
        str(run_dir.resolve()),
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


def build_snakemake_unlock_cmd(
    config_path: Path,
    snakefile: Path,
    workdir: Path | None = None,
) -> list[str]:
    run_dir = Path(workdir) if workdir is not None else config_path.parent
    return [
        sys.executable,
        "-m",
        "snakemake",
        "--snakefile",
        str(snakefile),
        "--directory",
        str(run_dir.resolve()),
        "--configfile",
        str(config_path),
        "--unlock",
    ]


def unlock_snakemake_run(config_path: Path, snakefile: Path, workdir: Path | None = None) -> int:
    cmd = build_snakemake_unlock_cmd(config_path, snakefile, workdir=workdir)
    result = subprocess.run(cmd, check=False)
    return result.returncode


def run_snakemake_target(
    config_path: Path,
    cores: int,
    target: str,
    snake_args: str,
    snakefile: Path,
    workdir: Path | None = None,
) -> int:
    cmd = build_snakemake_cmd(config_path, cores, target, snake_args, snakefile, workdir=workdir)
    result = subprocess.run(cmd, check=False)
    return result.returncode


def resume_command(outdir: Path) -> str:
    return f"babappasnake --resume --outdir {shlex.quote(str(Path(outdir).resolve()))}"


def resume_with_cds_command(outdir: Path) -> str:
    return f"{resume_command(outdir)} --cds /path/to/orthogroup_cds.fasta"


def print_guided_resume_banner(outdir: Path) -> None:
    print("Guided mode enabled.")
    print("You can stop after any step; completed outputs will be reused.")
    print(f"Resume command: {resume_command(outdir)}")


def print_cds_resume_instructions(outdir: Path) -> None:
    expected = (outdir / "user_supplied" / "orthogroup_cds.fasta").resolve()
    waiting_note = outdir / "orthogroup" / "WAITING_FOR_CDS.txt"
    print("")
    print("CDS checkpoint reached.")
    print(f"Place the orthogroup CDS FASTA at: {expected}")
    print(f"Then resume with: {resume_command(outdir)}")
    print("Or stage a CDS file directly with:")
    print(f"  {resume_with_cds_command(outdir)}")
    if waiting_note.exists():
        print(f"Detailed instructions were written to: {waiting_note}")


def staged_cds_ready(outdir: Path) -> bool:
    staged_cds = outdir / "user_supplied" / "orthogroup_cds.fasta"
    return staged_cds.exists() and staged_cds.stat().st_size > 0


def print_resume_after_interruption(outdir: Path) -> None:
    print("")
    print("The run did not complete.")
    print(f"After fixing the cause, resume from existing outputs with: {resume_command(outdir)}")


def build_step_plan(
    have_cds: bool,
    methods: list[str],
    trim_states: list[str],
    recombination_mode: str,
    run_asr: bool = True,
    tree_mode: str = "iqtree",
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

    steps.extend(
        [
            StepSpec(
                "tree_all_pathways",
                (
                    "Stage user-supplied tree for all selected pathways."
                    if parse_tree_mode(tree_mode) == "user"
                    else "Infer ML phylogeny with IQ-TREE for all selected pathways."
                ),
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
    )
    if run_asr:
        insert_at = next(i for i, step in enumerate(steps) if step.rule == "final_summary_all_pathways")
        steps[insert_at:insert_at] = [
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
        alias_at = next(i for i, step in enumerate(steps) if step.rule == "write_run_provenance")
        steps.insert(
            alias_at,
            StepSpec(
                "asr_primary_alias",
                "Write legacy top-level ASR completion alias.",
                ("asr/asr_done.json",),
            ),
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
    if all(existing):
        print("All expected outputs for this step already exist. Auto-skipping.")
        print_step_outputs(outdir, step.outputs)
        return 0
    can_skip = all(existing) or force_can_skip
    if can_skip:
        if force_can_skip:
            print("This step is skippable in the current context.")
    action = prompt_step_action(can_skip=can_skip)
    if action == "stop":
        print("Guided run stopped by user.")
        print(f"Resume command: {resume_command(outdir)}")
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
    summary = outdir / "orthogroup" / "orthogroup_summary.tsv"
    if not summary.exists():
        print("Orthogroup membership report unavailable (orthogroup_summary.tsv missing).")
        return

    included: list[tuple[str, str]] = []
    omitted: list[str] = []
    try:
        with open(summary, "r", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                species = str(row.get("species", "")).strip() or "(unknown)"
                selected = str(row.get("selected_members", "")).strip()
                if selected:
                    included.append((species, selected))
                else:
                    omitted.append(species)
    except Exception as exc:
        print(f"Orthogroup membership report unavailable (parse error: {exc}).")
        return

    print("\nOrthogroup membership summary:")
    print(f"  Included species with retained member(s): {len(included)}")
    for species, selected in included:
        print(f"    - {species}: {selected}")
    print(f"  Omitted species (no retained member in selected orthogroup/mode): {len(omitted)}")
    for species in omitted:
        print(f"    - {species}")


def prompt_for_cds_after_orthogroup(args: argparse.Namespace, outdir: Path) -> bool:
    user_cds = outdir / "user_supplied" / "orthogroup_cds.fasta"
    existing = user_cds.exists() and user_cds.stat().st_size > 0
    default = str(args.cds) if args.cds else (str(user_cds) if existing else "")
    print("")
    print("CDS input step (required for codon/tree/hyphy/codeml downstream).")
    print(f"Expected CDS checkpoint path: {user_cds.resolve()}")
    print(f"Resume command after placing CDS there: {resume_command(outdir)}")
    while True:
        entered = prompt_text(
            "Provide CDS FASTA path now (or press Enter to continue without CDS)",
            default,
            required=False,
        ).strip()
        if not entered:
            if staged_cds_ready(outdir):
                return True
            print("")
            print("No CDS file staged yet. BABAPPASnake will write the waiting note and stop at this checkpoint.")
            print(f"When the CDS is ready, place it at: {user_cds.resolve()}")
            print(f"Then resume with: {resume_command(outdir)}")
            return False

        candidate = Path(entered).expanduser()
        if not candidate.is_file():
            print(f"CDS file not found: {candidate}")
            default = ""
            continue
        if candidate.stat().st_size == 0:
            print(f"CDS file is empty: {candidate}")
            default = ""
            continue

        user_cds.parent.mkdir(parents=True, exist_ok=True)
        if candidate.resolve() != user_cds.resolve():
            shutil.copy2(candidate, user_cds)
        print(f"Staged CDS file: {user_cds}")
        print(f"Resume command for this run: {resume_command(outdir)}")
        return True


def set_outgroup_in_config(config_path: Path, outgroup_query: str) -> None:
    with open(config_path, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh) or {}
    cfg["outgroup_query"] = outgroup_query
    with open(config_path, "w", encoding="utf-8") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)


def set_tree_in_config(config_path: Path, tree_mode: str, user_tree: str | Path = "") -> None:
    with open(config_path, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh) or {}
    cfg["tree_mode"] = parse_tree_mode(tree_mode)
    cfg["user_tree"] = str(user_tree)
    with open(config_path, "w", encoding="utf-8") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)


def prompt_for_user_tree_after_cds(args: argparse.Namespace, outdir: Path) -> Path | None:
    staged_tree = outdir / "user_supplied" / "user_tree.nwk"
    existing = staged_tree.exists() and staged_tree.stat().st_size > 0
    default = str(args.user_tree) if getattr(args, "user_tree", None) else (str(staged_tree) if existing else "")
    while True:
        entered = prompt_text(
            "Provide Newick tree path now (or enter no to let IQ-TREE build trees)",
            default,
            required=False,
        ).strip()
        if not entered or entered.lower() in {"n", "no"}:
            return None

        candidate = Path(entered).expanduser()
        if not candidate.is_file():
            print(f"Tree file not found: {candidate}")
            default = ""
            continue

        staged_tree.parent.mkdir(parents=True, exist_ok=True)
        if candidate.resolve() != staged_tree.resolve():
            shutil.copy2(candidate, staged_tree)
        args.tree_mode = "user"
        args.user_tree = str(staged_tree)
        return staged_tree


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


def resume_state_path(outdir: Path) -> Path:
    return outdir / ".babappasnake" / "resume_state.json"


def load_resume_state(outdir: Path) -> dict:
    path = resume_state_path(outdir)
    if not path.exists() or path.stat().st_size == 0:
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def write_resume_state(outdir: Path, state: dict) -> None:
    path = resume_state_path(outdir)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")


def validate_resume_request(args: argparse.Namespace, argv: list[str]) -> None:
    allowed = {
        "--resume",
        "--outdir",
        "--cds",
        "--outgroup",
        "--tree",
        "--user-tree",
        "--tree-mode",
        "--threads",
        "--guided",
        "--snake-args",
        "--interactive",
    }
    disallowed = {
        "--prot",
        "--query",
        "--orthogroup-method",
        "--orthogroup-proteins",
        "--orthology-mode",
        "--coverage",
        "--alignment-methods",
        "--trim-strategy",
        "--use-clipkit",
        "--clipkit-mode-protein",
        "--clipkit-mode-codon",
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
        "--absrel-p",
        "--absrel-dynamic-start",
        "--absrel-dynamic-step",
        "--absrel-dynamic-max",
        "--meme-p",
    }
    seen = {token for token in argv if token.startswith("--")}
    blocked = sorted(seen & disallowed)
    if blocked:
        raise SystemExit(
            "Resume can only override outdir, cds, outgroup, tree, tree-mode, threads, guided, and snake-args. "
            "Start a new run directory to change analysis settings: " + ", ".join(blocked)
        )
    unknown = sorted(token for token in seen if token not in allowed and token not in disallowed)
    if unknown:
        return
    if not getattr(args, "outdir", ""):
        raise SystemExit("--resume requires --outdir")


def _argv_has(argv: list[str], flag: str) -> bool:
    return flag in argv


def maybe_prompt_outgroup_after_cds(args: argparse.Namespace, config_path: Path, have_cds: bool) -> None:
    if not have_cds:
        return
    if not is_tty_interactive():
        return
    outgroup_prompt = prompt_text(
        "Optional outgroup query for rooting (press Enter to skip)",
        str(getattr(args, "outgroup", "") or ""),
        required=False,
    ).strip()
    args.outgroup = outgroup_prompt
    set_outgroup_in_config(config_path, args.outgroup)


def prepare_resume_run(args: argparse.Namespace, argv: list[str]) -> tuple[argparse.Namespace, Path, Path]:
    validate_resume_request(args, argv)
    outdir = Path(args.outdir)
    config_path = outdir / "config.yaml"
    if not config_path.exists():
        raise SystemExit(f"Cannot resume because config.yaml is missing: {config_path}")

    requested_user_tree = getattr(args, "user_tree", None)
    requested_tree_mode = getattr(args, "tree_mode", "iqtree")

    with open(config_path, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh) or {}

    methods = [str(x) for x in cfg.get("alignment_methods", ["babappalign", "mafft", "prank"])]
    trim_states = [str(x) for x in cfg.get("trim_states", ["raw", "clipkit"])]

    args.prot = cfg.get("proteomes_dir", getattr(args, "prot", None))
    args.query = cfg.get("query_fasta", getattr(args, "query", None))
    args.orthogroup_proteins = cfg.get("external_orthogroup_proteins", getattr(args, "orthogroup_proteins", None))
    args.orthogroup_method = "orthofinder"
    args.orthology_mode = str(cfg.get("orthology_mode", "representative"))
    args.alignment_methods = str(cfg.get("alignment_method_option", alignment_option_from_methods(methods)))
    args.coverage = float(cfg.get("coverage", getattr(args, "coverage", 0.70)))
    args.trim_strategy = str(cfg.get("trim_strategy", "both"))
    args.use_clipkit = "yes" if bool(cfg.get("use_clipkit", True)) else "no"
    args.clipkit_mode_protein = str(cfg.get("clipkit_mode_protein", "kpic-smart-gap"))
    args.clipkit_mode_codon = str(cfg.get("clipkit_mode_codon", "kpic-smart-gap"))
    args.iqtree_bootstrap = int(cfg.get("iqtree_bootstrap", 1000))
    args.iqtree_bnni = "yes" if bool(cfg.get("iqtree_bnni", False)) else "no"
    args.iqtree_model = str(cfg.get("iqtree_model", "MFP"))
    args.absrel_branches = str(cfg.get("absrel_branches", "Leaves"))
    args.meme_branches = str(cfg.get("meme_branches", "Leaves"))
    args.codeml_codonfreq = int(cfg.get("codeml_codonfreq", 7))
    args.run_asr = "yes" if bool(cfg.get("run_asr", True)) else "no"
    args.recombination = str(cfg.get("recombination", "none"))
    args.gard_mode = str(cfg.get("gard_mode", "Faster"))
    args.gard_rate_classes = int(cfg.get("gard_rate_classes", 3))
    args.absrel_p = float(cfg.get("absrel_p", 0.05))
    args.absrel_dynamic_start = float(cfg.get("absrel_dynamic_start", 0.05))
    args.absrel_dynamic_step = float(cfg.get("absrel_dynamic_step", 0.01))
    args.absrel_dynamic_max = float(cfg.get("absrel_dynamic_max", 0.2))
    args.meme_p = float(cfg.get("meme_p", 0.05))

    if _argv_has(argv, "--threads"):
        cfg["threads"] = int(args.threads)
    else:
        args.threads = int(cfg.get("threads", getattr(args, "threads", DEFAULT_TOTAL_THREADS)))
    if _argv_has(argv, "--guided"):
        cfg["guided_mode"] = str(args.guided).strip().lower() == "yes"
    else:
        args.guided = "yes" if bool(cfg.get("guided_mode", False)) else "no"
    if _argv_has(argv, "--snake-args"):
        cfg["snake_args"] = str(args.snake_args or "")
    else:
        args.snake_args = str(cfg.get("snake_args", ""))
    if _argv_has(argv, "--outgroup"):
        cfg["outgroup_query"] = str(args.outgroup or "").strip()
    else:
        args.outgroup = str(cfg.get("outgroup_query", ""))
    args.tree_mode = str(cfg.get("tree_mode", "iqtree"))
    args.user_tree = str(cfg.get("user_tree", ""))
    tree_override = _argv_has(argv, "--tree") or _argv_has(argv, "--user-tree")
    tree_mode_override = _argv_has(argv, "--tree-mode")
    if tree_override:
        candidate = Path(str(requested_user_tree)).expanduser()
        if not candidate.is_file():
            raise SystemExit(f"--tree must point to an existing Newick tree file: {candidate}")
        user_tree_dst = outdir / "user_supplied" / "user_tree.nwk"
        user_tree_dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(candidate, user_tree_dst)
        cfg["tree_mode"] = "user"
        cfg["user_tree"] = str(user_tree_dst.resolve())
        args.tree_mode = "user"
        args.user_tree = str(user_tree_dst)
    elif tree_mode_override:
        requested_tree_mode_normalized = parse_tree_mode(requested_tree_mode)
        if requested_tree_mode_normalized == "user" and not str(cfg.get("user_tree", "")).strip():
            raise SystemExit("--tree-mode user on resume requires --tree unless a user tree is already staged.")
        cfg["tree_mode"] = requested_tree_mode_normalized
        if requested_tree_mode_normalized == "iqtree":
            cfg["user_tree"] = ""
        args.tree_mode = requested_tree_mode_normalized
        args.user_tree = str(cfg.get("user_tree", ""))
    else:
        args.tree_mode = str(cfg.get("tree_mode", "iqtree"))
        args.user_tree = str(cfg.get("user_tree", ""))

    user_cds = outdir / "user_supplied" / "orthogroup_cds.fasta"
    if _argv_has(argv, "--cds") and args.cds:
        candidate = Path(args.cds)
        if not candidate.is_file():
            raise SystemExit(f"--cds must point to an existing FASTA file: {candidate}")
        user_cds.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(candidate, user_cds)
    cfg["user_cds"] = str(user_cds.resolve())

    method_count = max(1, len(methods))
    pathway_count = max(1, len(methods) * len(trim_states))
    cfg["orthogroup_method"] = "orthofinder"
    cfg["orthogroup_source"] = str(
        cfg.get("orthogroup_source")
        or ("external" if cfg.get("external_orthogroup_proteins") else "orthofinder")
    )
    cfg["orthology_mode"] = args.orthology_mode
    cfg["threads"] = int(args.threads)
    cfg["per_method_cores"] = compute_per_method_cores(int(args.threads), method_count)
    cfg["per_pathway_cores"] = compute_per_pathway_cores(int(args.threads), pathway_count)
    cfg["guided_mode"] = str(args.guided).strip().lower() == "yes"
    cfg["snake_args"] = str(args.snake_args or "")
    cfg["outgroup_query"] = str(args.outgroup or "").strip()
    cfg["tree_mode"] = parse_tree_mode(str(cfg.get("tree_mode", "iqtree")))
    cfg["user_tree"] = str(cfg.get("user_tree", ""))

    with open(config_path, "w", encoding="utf-8") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)
    write_resume_state(
        outdir,
        {
            "guided_mode": cfg["guided_mode"],
            "snake_args": cfg["snake_args"],
            "threads": cfg["threads"],
            "orthogroup_method": "orthofinder",
            "orthogroup_source": cfg["orthogroup_source"],
            "orthology_mode": args.orthology_mode,
            "tree_mode": cfg["tree_mode"],
            "user_tree": cfg["user_tree"],
        },
    )
    return args, config_path, outdir


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
    source = orthogroup_source_from_args(args)
    step_description = (
        "Stage externally curated orthogroup proteins."
        if source == "external"
        else f"Build orthogroup proteins with OrthoFinder ({args.orthology_mode} mode)."
    )
    orthogroup_step = StepSpec(
        "define_orthogroup",
        step_description,
        (
            "orthogroup/orthogroup_proteins.fasta",
            "orthogroup/orthogroup_headers.txt",
            "orthogroup/orthogroup_summary.tsv",
            "orthogroup/orthogroup_metadata.json",
        ),
    )
    print_guided_resume_banner(outdir)
    code = run_guided_step(
        config_path=config_path,
        outdir=outdir,
        cores=cores,
        snake_args=snake_args,
        snakefile=snakefile,
        step=orthogroup_step,
        index=1,
        total=None,
    )
    if code != 0:
        return code

    print_orthogroup_membership(outdir)

    have_cds = staged_cds_ready(outdir)
    if is_tty_interactive():
        have_cds = prompt_for_cds_after_orthogroup(args, outdir)
        if have_cds:
            staged_tree = outdir / "user_supplied" / "user_tree.nwk"
            if tree_mode_from_args(args) == "user" and staged_tree.exists() and staged_tree.stat().st_size > 0:
                args.tree_mode = "user"
                args.user_tree = str(staged_tree)
                set_tree_in_config(config_path, "user", staged_tree)
                print(f"Using user-supplied tree: {staged_tree}")
            elif getattr(args, "tree_choice_confirmed", False):
                args.tree_mode = "iqtree"
                args.user_tree = ""
                set_tree_in_config(config_path, "iqtree", "")
                print("Using IQ-TREE for pathway tree inference based on the earlier interactive selection.")
            else:
                use_user_tree = prompt_yes_no(
                    "Use a user-supplied tree instead of building one with IQ-TREE?",
                    False,
                )
                if use_user_tree:
                    staged = prompt_for_user_tree_after_cds(args, outdir)
                    if staged:
                        set_tree_in_config(config_path, "user", staged)
                        print(f"Using user-supplied tree: {staged}")
                    else:
                        args.tree_mode = "iqtree"
                        args.user_tree = ""
                        set_tree_in_config(config_path, "iqtree", "")
                        print("No user tree supplied; IQ-TREE will infer pathway trees.")
                else:
                    args.tree_mode = "iqtree"
                    args.user_tree = ""
                    set_tree_in_config(config_path, "iqtree", "")
                    print("No user tree supplied; IQ-TREE will infer pathway trees.")
        else:
            args.outgroup = str(args.outgroup or "").strip()
            set_outgroup_in_config(config_path, args.outgroup)
        if have_cds:
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

    steps = build_step_plan(
        have_cds,
        methods,
        trim_states,
        args.recombination,
        run_asr=parse_yes_no_bool(getattr(args, "run_asr", "yes"), True),
        tree_mode=tree_mode_from_args(args),
    )
    total_steps = 1 + len(steps)
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
    args = parse_args()
    if args.resume:
        args, config_path, outdir = prepare_resume_run(args, sys.argv[1:])
        snakefile = ir.files("babappasnake").joinpath("workflow", "Snakefile")
        methods = parse_alignment_method_option(args.alignment_methods)
        trim_states = trim_states_from_strategy(args.trim_strategy)
        have_cds = staged_cds_ready(outdir)
        unlock_rc = unlock_snakemake_run(config_path, snakefile, workdir=outdir)
        if unlock_rc != 0:
            print(
                "[WARN] Snakemake unlock returned a non-zero status; continuing resume attempt.",
                file=sys.stderr,
            )
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
        if rc == 0 and not have_cds and not staged_cds_ready(outdir):
            print_cds_resume_instructions(outdir)
        if rc != 0:
            print_resume_after_interruption(outdir)
        raise SystemExit(rc)

    args = maybe_prompt_interactive(args)
    args.tree_mode = tree_mode_from_args(args)
    validate_inputs(args)
    args.orthogroup_method = "orthofinder"
    args.orthology_mode = str(args.orthology_mode).strip().lower()
    orthogroup_source = orthogroup_source_from_args(args)
    args.recombination = parse_recombination_mode(args.recombination)
    recombination_mode = effective_recombination_mode(args.recombination)
    tree_mode = parse_tree_mode(args.tree_mode)
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
    print(f"[INFO] Orthogroup source: {orthogroup_source}.")
    if orthogroup_source == "orthofinder":
        print(f"[INFO] Orthology mode: {args.orthology_mode}.")
    print(f"[INFO] Tree mode: {tree_mode}.")
    if tree_mode == "user":
        print(f"[INFO] User-supplied tree: {args.user_tree}.")
    if recombination_mode == "gard":
        print(
            f"[INFO] GARD parameters: mode={args.gard_mode}, "
            f"rate_classes={args.gard_rate_classes}."
        )

    resolved, missing = resolve_tools()
    missing = filter_missing_tools_for_run(missing, tree_mode)
    if missing:
        print(format_missing_tools(missing), file=sys.stderr)
        raise SystemExit(2)
    if orthogroup_source == "orthofinder":
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
    snakefile = ir.files("babappasnake").joinpath("workflow", "Snakefile")

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
    if rc == 0 and not staged_cds_ready(outdir):
        print_cds_resume_instructions(outdir)
    if rc != 0:
        print_resume_after_interruption(outdir)
    raise SystemExit(rc)


if __name__ == "__main__":
    main()
