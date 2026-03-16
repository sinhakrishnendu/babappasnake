from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from importlib import resources
from pathlib import Path
from typing import Any

import yaml

from ._version import __version__


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="babappasnake",
        description=(
            "Run the packaged BABAPPAsnake Snakemake workflow. "
            "The direct run form is: babappasnake --input PROTEOMES --query QUERY --outgroup TAXON "
            "--clipkit yes|no --iqtree protein|cds"
        ),
    )
    parser.add_argument("--input", help="Path to the proteome directory.")
    parser.add_argument("--query", help="Path to the query FASTA.")
    parser.add_argument("--outgroup", help="Exact outgroup taxon name required in the final orthogroup.")
    parser.add_argument(
        "--clipkit",
        choices=["yes", "no"],
        help="Enable or disable ClipKIT trimming.",
    )
    parser.add_argument(
        "--iqtree",
        choices=["protein", "cds"],
        help="Choose whether IQ-TREE runs on the protein or CDS alignment.",
    )
    parser.add_argument("--config", help="Optional base YAML config file to merge with CLI overrides.")
    parser.add_argument(
        "--threads",
        type=int,
        default=max(os.cpu_count() or 1, 1),
        help="Snakemake core count. Default: detected CPU count.",
    )
    parser.add_argument(
        "--output",
        default="results",
        help="Results directory written into the runtime config. Default: results",
    )
    parser.add_argument(
        "--coverage-thresholds",
        help="Comma-separated RBH coverage thresholds, for example 50,60,70,80,90,95.",
    )
    parser.add_argument(
        "--cds-input",
        help="Optional CDS FASTA to seed the checkpoint request automatically.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume an existing run. Snakemake resumes incomplete work by design; this flag is informational.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Render the Snakemake DAG and commands without executing tools.",
    )
    parser.add_argument(
        "--use-conda",
        action="store_true",
        help="Forward --use-conda to Snakemake.",
    )
    parser.add_argument(
        "--conda-prefix",
        help="Optional conda environment prefix forwarded to Snakemake.",
    )
    parser.add_argument(
        "--init-config",
        metavar="PATH",
        help="Write the bundled default config to PATH and exit.",
    )
    parser.add_argument(
        "--validate-installation",
        action="store_true",
        help="Validate that the package can locate the bundled workflow assets and Snakemake.",
    )
    parser.add_argument(
        "--summarize",
        metavar="RESULTS_DIR",
        help="Print RESULTS_DIR/reports/executive_summary.txt and exit.",
    )
    parser.add_argument(
        "--snakemake-arg",
        action="append",
        default=[],
        help="Additional argument forwarded to Snakemake. Repeat this flag for multiple arguments.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return parser


def repo_workflow_root() -> Path | None:
    for parent in Path(__file__).resolve().parents:
        if (parent / "Snakefile").exists() and (parent / "workflow").is_dir() and (parent / "envs").is_dir():
            return parent
    return None


def copy_traversable_tree(source: Any, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    for child in source.iterdir():
        child_destination = destination / child.name
        if child.is_dir():
            copy_traversable_tree(child, child_destination)
        else:
            child_destination.parent.mkdir(parents=True, exist_ok=True)
            with resources.as_file(child) as resolved_child:
                shutil.copy2(resolved_child, child_destination)


def bundled_workflow_root() -> Path:
    cache_candidates = []
    xdg_cache_home = os.environ.get("XDG_CACHE_HOME", "").strip()
    if xdg_cache_home:
        cache_candidates.append(Path(xdg_cache_home) / "babappasnake" / "workflow_bundle" / __version__)
    cache_candidates.append(Path.home() / ".cache" / "babappasnake" / "workflow_bundle" / __version__)
    cache_candidates.append(Path(tempfile.gettempdir()) / "babappasnake" / "workflow_bundle" / __version__)

    bundle_root = resources.files("babappasnake").joinpath("_bundle")
    last_error = None
    for cache_root in cache_candidates:
        try:
            marker = cache_root / ".bundle_ready"
            if marker.exists():
                return cache_root

            if cache_root.exists():
                shutil.rmtree(cache_root)
            cache_root.mkdir(parents=True, exist_ok=True)
            copy_traversable_tree(bundle_root, cache_root)
            marker.write_text(__version__ + "\n")
            return cache_root
        except PermissionError as exc:
            last_error = exc
            continue

    if last_error is not None:
        raise last_error
    raise RuntimeError("Unable to materialize bundled workflow assets.")


def workflow_root() -> Path:
    source_root = repo_workflow_root()
    if source_root is not None:
        return source_root
    return bundled_workflow_root()


def load_yaml_config(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text())
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ValueError(f"Config file must contain a mapping at the top level: {path}")
    return data


def ensure_nested_dict(root: dict[str, Any], key: str) -> dict[str, Any]:
    value = root.get(key)
    if not isinstance(value, dict):
        value = {}
        root[key] = value
    return value


def parse_coverage_thresholds(raw: str | None) -> list[int] | None:
    if not raw:
        return None
    values = []
    for item in raw.split(","):
        stripped = item.strip()
        if not stripped:
            continue
        values.append(int(stripped))
    if not values:
        raise ValueError("coverage thresholds string did not contain any integers.")
    if sorted(values) != values:
        raise ValueError("coverage thresholds must be sorted in ascending order.")
    return values


def resolve_relative_config_path(config_dir: Path, value: Any) -> Any:
    if not isinstance(value, str) or not value.strip():
        return value
    candidate = Path(value).expanduser()
    if candidate.is_absolute():
        return str(candidate.resolve())
    resolved = (config_dir / candidate).resolve()
    if resolved.exists():
        return str(resolved)
    return value


def create_runtime_config(args: argparse.Namespace, root: Path) -> tuple[dict[str, Any], Path]:
    base_config_path = Path(args.config).expanduser().resolve() if args.config else root / "config" / "config.yaml"
    config = load_yaml_config(base_config_path)
    config_dir = base_config_path.parent
    for path_key in ("species_metadata_tsv", "cds_input_fasta"):
        if path_key in config:
            config[path_key] = resolve_relative_config_path(config_dir, config[path_key])
    config["proteomes_dir"] = str(Path(args.input).expanduser().resolve())
    config["query_fasta"] = str(Path(args.query).expanduser().resolve())
    config["outgroup_name"] = args.outgroup
    config["results_dir"] = str(Path(args.output).expanduser().resolve())
    if args.cds_input:
        config["cds_input_fasta"] = str(Path(args.cds_input).expanduser().resolve())

    clipkit_config = ensure_nested_dict(config, "clipkit")
    clipkit_config["enabled"] = args.clipkit == "yes"

    iqtree_config = ensure_nested_dict(config, "iqtree")
    iqtree_config["input_type"] = args.iqtree

    alignment_config = ensure_nested_dict(config, "alignment")
    if not alignment_config.get("model_path"):
        default_model = (Path.home() / ".cache" / "babappalign" / "models" / "babappascore.pt").resolve()
        if default_model.exists():
            alignment_config["model_path"] = str(default_model)

    thresholds = parse_coverage_thresholds(args.coverage_thresholds)
    if thresholds is not None:
        config["coverage_thresholds"] = thresholds

    runtime_dir = Path(config["results_dir"]) / ".babappasnake"
    runtime_dir.mkdir(parents=True, exist_ok=True)
    runtime_config_path = runtime_dir / "runtime_config.yaml"
    runtime_config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    return config, runtime_config_path


def snakemake_command(args: argparse.Namespace, root: Path, runtime_config_path: Path) -> list[str]:
    command = [
        sys.executable,
        "-m",
        "snakemake",
        "--snakefile",
        str(root / "Snakefile"),
        "--directory",
        str(Path.cwd()),
        "--configfile",
        str(runtime_config_path),
        "--cores",
        str(args.threads),
        "--rerun-incomplete",
    ]
    if args.use_conda:
        command.append("--use-conda")
    if args.conda_prefix:
        command.extend(["--conda-prefix", args.conda_prefix])
    if args.dry_run:
        command.append("--dry-run")
    for extra_arg in args.snakemake_arg:
        command.append(extra_arg)
    return command


def validate_run_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    required_flags = {
        "--input": args.input,
        "--query": args.query,
        "--outgroup": args.outgroup,
        "--clipkit": args.clipkit,
        "--iqtree": args.iqtree,
    }
    missing = [flag for flag, value in required_flags.items() if not value]
    if missing:
        parser.error(
            "the following arguments are required for workflow execution: " + ", ".join(missing)
        )
    if args.threads < 1:
        parser.error("--threads must be at least 1")


def handle_init_config(path_text: str) -> int:
    root = workflow_root()
    output_path = Path(path_text).expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(root / "config" / "config.yaml", output_path)
    print(output_path)
    return 0


def handle_validate_installation() -> int:
    root = workflow_root()
    snakefile = root / "Snakefile"
    if not snakefile.exists():
        raise FileNotFoundError(f"Bundled Snakefile was not found at {snakefile}")
    import snakemake  # noqa: F401

    print(f"babappasnake {__version__}")
    print(f"workflow_root={root}")
    print(f"snakefile={snakefile}")
    return 0


def handle_summarize(results_dir: str) -> int:
    summary_path = Path(results_dir).expanduser().resolve() / "reports" / "executive_summary.txt"
    if not summary_path.exists():
        raise FileNotFoundError(f"Executive summary not found: {summary_path}")
    print(summary_path.read_text(), end="")
    return 0


def run_workflow(parser: argparse.ArgumentParser, args: argparse.Namespace) -> int:
    validate_run_args(parser, args)
    root = workflow_root()
    config, runtime_config_path = create_runtime_config(args, root)
    command = snakemake_command(args, root, runtime_config_path)
    runtime_dir = Path(config["results_dir"]) / ".babappasnake"
    cache_dir = runtime_dir / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["XDG_CACHE_HOME"] = str(cache_dir)
    env["HOME"] = str(runtime_dir)
    print(" ".join(command))
    return subprocess.run(command, env=env).returncode


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.init_config:
        return handle_init_config(args.init_config)
    if args.validate_installation:
        return handle_validate_installation()
    if args.summarize:
        return handle_summarize(args.summarize)
    return run_workflow(parser, args)
