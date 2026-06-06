from __future__ import annotations

import shutil
from dataclasses import dataclass


@dataclass(frozen=True)
class ToolSpec:
    key: str
    choices: tuple[str, ...]
    required: bool = True
    install_hint: str = ""


TOOL_SPECS: tuple[ToolSpec, ...] = (
    ToolSpec("python", ("python",), True, "Provided by the active environment."),
    ToolSpec("blastp", ("blastp",), False, "conda install -c conda-forge -c bioconda blast"),
    ToolSpec("makeblastdb", ("makeblastdb",), False, "conda install -c conda-forge -c bioconda blast"),
    ToolSpec("orthofinder", ("orthofinder",), False, "conda install -c conda-forge -c bioconda orthofinder"),
    ToolSpec("iqtree", ("iqtree2", "iqtree3", "iqtree"), True, "conda install -c conda-forge -c bioconda iqtree"),
    ToolSpec("hyphy", ("hyphy",), True, "conda install -c conda-forge -c bioconda hyphy"),
    ToolSpec("codeml", ("codeml",), True, "conda install -c conda-forge -c bioconda paml"),
    ToolSpec("clipkit", ("clipkit",), False, "conda install -c conda-forge -c bioconda clipkit"),
    ToolSpec("babappalign", ("babappalign",), False, "pip install babappalign"),
    ToolSpec("mafft", ("mafft",), False, "conda install -c conda-forge -c bioconda mafft"),
    ToolSpec("prank", ("prank",), False, "conda install -c conda-forge -c bioconda prank"),
    ToolSpec("pal2nal", ("pal2nal.pl", "pal2nal"), False, "conda install -c conda-forge -c bioconda pal2nal"),
)


def resolve_tools() -> tuple[dict[str, str], list[tuple[ToolSpec, tuple[str, ...]]]]:
    resolved: dict[str, str] = {}
    missing: list[tuple[ToolSpec, tuple[str, ...]]] = []
    for spec in TOOL_SPECS:
        found = None
        for exe in spec.choices:
            path = shutil.which(exe)
            if path:
                found = path
                break
        if found:
            resolved[spec.key] = found
        elif spec.required:
            missing.append((spec, spec.choices))
    return resolved, missing


def format_missing_tools(missing: list[tuple[ToolSpec, tuple[str, ...]]]) -> str:
    lines = [
        "Missing required external tools. Install them before running babappasnake:",
        "",
    ]
    for spec, choices in missing:
        lines.append(f"- {spec.key}: expected one of {', '.join(choices)}")
        if spec.install_hint:
            lines.append(f"  install: {spec.install_hint}")
    lines.extend([
        "",
        "Tip: on Apple Silicon, IQ-TREE may be installed as iqtree3; babappasnake accepts iqtree2, iqtree3, or iqtree automatically.",
    ])
    return "\n".join(lines)
