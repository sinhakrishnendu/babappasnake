import csv
import json
import os
import re
import shutil
from pathlib import Path

from snakemake.utils import min_version

configfile: "config/config.yaml"

min_version("7.0")
shell.executable("/bin/bash")

PROJECT_ROOT = Path(workflow.basedir).resolve()
WORKFLOW_DIR = (PROJECT_ROOT / "workflow").as_posix()
SCRIPT_DIR = (PROJECT_ROOT / "workflow" / "scripts").as_posix()
ENV_DIR = (PROJECT_ROOT / "envs").as_posix()

RESULTS_DIR = Path(config.get("results_dir", "results")).as_posix()
LOGS_DIR = f"{RESULTS_DIR}/logs"
VALIDATION_DIR = f"{RESULTS_DIR}/validation"
RBH_DIR = f"{RESULTS_DIR}/rbh"
THRESHOLD_DIR = f"{RESULTS_DIR}/threshold_comparison"
SELECTED_DIR = f"{RESULTS_DIR}/selected_orthogroup"
CHECKPOINT_DIR = f"{RESULTS_DIR}/cds_input_checkpoint"
ALIGNMENT_DIR = f"{RESULTS_DIR}/alignment"
ALIGNMENT_BAB_DIR = f"{ALIGNMENT_DIR}/babappalign"
ALIGNMENT_BAB_PROTEIN_DIR = f"{ALIGNMENT_BAB_DIR}/protein"
ALIGNMENT_BAB_CDS_DIR = f"{ALIGNMENT_BAB_DIR}/cds"
ALIGNMENT_CLIPKIT_DIR = f"{ALIGNMENT_DIR}/clipkit"
ALIGNMENT_CLIPKIT_PROTEIN_DIR = f"{ALIGNMENT_CLIPKIT_DIR}/protein"
ALIGNMENT_CLIPKIT_CDS_DIR = f"{ALIGNMENT_CLIPKIT_DIR}/cds"
IQTREE_DIR = f"{RESULTS_DIR}/iqtree"
HYPHY_DIR = f"{RESULTS_DIR}/hyphy"
CODEML_DIR = f"{RESULTS_DIR}/codeml"
REPORT_DIR = f"{RESULTS_DIR}/reports"

PROTEOME_EXTENSIONS = tuple(
    suffix.lower()
    for suffix in config.get(
        "proteome_extensions",
        [".faa", ".fa", ".fasta", ".fsa", ".pep", ".protein"],
    )
)
THRESHOLDS = [int(value) for value in config["coverage_thresholds"]]
if THRESHOLDS != sorted(THRESHOLDS):
    raise ValueError("coverage_thresholds must be sorted in ascending order.")


def as_bool_config(value):
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def require_choice(name, value, allowed):
    normalized = str(value).strip().lower()
    if normalized not in allowed:
        raise ValueError(f"{name} must be one of {sorted(allowed)}, got {value!r}")
    return normalized


def sanitize_name(value):
    sanitized = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")
    return sanitized or "item"


def discover_proteomes():
    proteome_dir = Path(config["proteomes_dir"]).expanduser()
    if not proteome_dir.exists():
        raise ValueError(f"Proteome directory does not exist: {proteome_dir}")

    seen = {}
    records = []
    for path in sorted(proteome_dir.iterdir()):
        if not path.is_file():
            continue
        if path.suffix.lower() not in PROTEOME_EXTENSIONS:
            continue
        species_id = sanitize_name(path.stem)
        if species_id in seen:
            raise ValueError(
                f"Duplicate sanitized proteome identifier {species_id!r} "
                f"for {seen[species_id]!r} and {str(path)!r}."
            )
        seen[species_id] = str(path)
        records.append(
            {
                "species_id": species_id,
                "file_label": path.stem,
                "path": str(path),
            }
        )

    if not records:
        raise ValueError(
            f"No proteome FASTA files matching {PROTEOME_EXTENSIONS} were found in {proteome_dir}"
        )
    return records


PROTEOME_RECORDS = discover_proteomes()
SPECIES_IDS = [record["species_id"] for record in PROTEOME_RECORDS]
PROTEOME_BY_ID = {record["species_id"]: record for record in PROTEOME_RECORDS}
PROTEOME_JSON = json.dumps(PROTEOME_RECORDS)

QUERY_FASTA = str(Path(config["query_fasta"]).expanduser())

FALLBACK_EXECUTABLES = {
    "blastp": [],
    "makeblastdb": [],
    "diamond": [],
    "clipkit": [],
    "babappalign": [],
    "hyphy": [],
    "iqtree3": [],
    "codeml": [],
}


def resolve_executable(tool_name):
    configured = config.get("executables", {}).get(tool_name, "")
    candidates = []
    if configured:
        candidates.append(os.path.expanduser(configured))
    candidates.extend(FALLBACK_EXECUTABLES.get(tool_name, []))
    candidates.append(tool_name)

    for candidate in candidates:
        if not candidate:
            continue
        expanded = os.path.expanduser(candidate)
        if shutil.which(expanded) or Path(expanded).exists():
            return expanded
    return candidates[-1]


def resolve_babappalign_model():
    configured = config.get("alignment", {}).get("model_path", "")
    candidates = []
    if configured:
        candidates.append(os.path.expanduser(configured))
    candidates.append(os.path.expanduser("~/.cache/babappalign/models/babappascore.pt"))
    for candidate in candidates:
        if candidate and Path(candidate).exists():
            return candidate
    return candidates[0] if configured else ""


TOOL_PATHS = {name: resolve_executable(name) for name in FALLBACK_EXECUTABLES}
BABAPPALIGN_MODEL = resolve_babappalign_model()
SPECIES_METADATA_TSV = config.get("species_metadata_tsv", "")
CLIPKIT_ENABLED = as_bool_config(config.get("clipkit", {}).get("enabled", True))
IQTREE_INPUT_TYPE = require_choice(
    "iqtree.input_type",
    config.get("iqtree", {}).get("input_type", "protein"),
    {"protein", "cds"},
)


def babappalign_protein_alignment():
    return f"{ALIGNMENT_BAB_PROTEIN_DIR}/aligned_proteins.fasta"


def babappalign_cds_alignment():
    return f"{ALIGNMENT_BAB_CDS_DIR}/aligned_cds.fasta"


def clipkit_trimmed_protein_alignment():
    return f"{ALIGNMENT_CLIPKIT_PROTEIN_DIR}/trimmed_proteins.fasta"


def clipkit_trimmed_cds_alignment():
    return f"{ALIGNMENT_CLIPKIT_CDS_DIR}/trimmed_cds.fasta"


def projected_clipkit_trimmed_cds_alignment():
    return f"{ALIGNMENT_CLIPKIT_CDS_DIR}/projected_trimmed_cds.fasta"


def effective_codon_alignment_path():
    if CLIPKIT_ENABLED:
        return clipkit_trimmed_cds_alignment()
    return babappalign_cds_alignment()


def effective_iqtree_alignment_path():
    if IQTREE_INPUT_TYPE == "protein":
        if CLIPKIT_ENABLED:
            return clipkit_trimmed_protein_alignment()
        return babappalign_protein_alignment()
    if CLIPKIT_ENABLED:
        return clipkit_trimmed_cds_alignment()
    return babappalign_cds_alignment()


def effective_iqtree_sequence_type():
    iqtree_cfg = config.get("iqtree", {})
    if IQTREE_INPUT_TYPE == "protein":
        return iqtree_cfg.get("protein_sequence_type") or iqtree_cfg.get("sequence_type", "AA")
    return iqtree_cfg.get("cds_sequence_type") or "DNA"


def effective_iqtree_model():
    iqtree_cfg = config.get("iqtree", {})
    if IQTREE_INPUT_TYPE == "protein":
        return iqtree_cfg.get("protein_model") or iqtree_cfg.get("model", "MFP")
    return iqtree_cfg.get("cds_model") or iqtree_cfg.get("model", "MFP")


def alignment_validation_dependencies(_wildcards):
    dependencies = [
        babappalign_protein_alignment(),
        babappalign_cds_alignment(),
        f"{SELECTED_DIR}/orthogroup_members.tsv",
    ]
    if CLIPKIT_ENABLED:
        dependencies.extend(
            [
                clipkit_trimmed_protein_alignment(),
                clipkit_trimmed_cds_alignment(),
                projected_clipkit_trimmed_cds_alignment(),
                f"{ALIGNMENT_CLIPKIT_CDS_DIR}/projection_summary.json",
            ]
        )
    return dependencies


def threshold_output_dir(threshold):
    return f"{THRESHOLD_DIR}/coverage_{int(threshold)}"


def threshold_member_file(threshold):
    return f"{threshold_output_dir(threshold)}/orthogroup_members.tsv"


def threshold_fasta_file(threshold):
    return f"{threshold_output_dir(threshold)}/orthogroup_proteins.faa"


def cds_checkpoint_dir(_wildcards):
    return checkpoints.cds_input_checkpoint.get().output[0]


def selected_codeml_branch_ids():
    selection_ckpt = checkpoints.select_codeml_targets.get()
    branch_ids = []
    with open(selection_ckpt.output.tsv) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        branch_ids = [row["branch_id"] for row in reader]
    return branch_ids


def selected_codeml_markers(_wildcards):
    return expand(
        f"{CODEML_DIR}/{{branch}}/branch_site.done",
        branch=selected_codeml_branch_ids(),
    )


rule all:
    input:
        f"{REPORT_DIR}/executive_summary.txt"


include: "workflow/rules/rbh.smk"
include: "workflow/rules/checkpoint.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/tree.smk"
include: "workflow/rules/hyphy.smk"
include: "workflow/rules/codeml.smk"
include: "workflow/rules/reporting.smk"
