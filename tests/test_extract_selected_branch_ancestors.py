from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO


METHOD = "babappalign"
TRIM = "clipkit"


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _write_fasta(path: Path, records: dict[str, str]) -> None:
    lines = []
    for k, v in records.items():
        lines.append(f">{k}")
        lines.append(v)
    _write(path, "\n".join(lines) + "\n")


def _seed_common(tmp_path: Path, tree_newick: str, branch_rows: list[tuple[str, str]]) -> Path:
    outdir = tmp_path / "run"

    _write(
        outdir / "branchsite" / METHOD / TRIM / "branchsite_results.tsv",
        "foreground_branch\tLRT\tp_value\tq_value\tsignificant_BH_0.05\n"
        + "".join(f"{fg}\t10.0\t0.001\t0.01\t{sig}\n" for fg, sig in branch_rows),
    )
    _write(outdir / "tree" / METHOD / TRIM / "orthogroup.rooted.treefile", tree_newick + "\n")
    _write_fasta(
        outdir / "alignments" / METHOD / TRIM / "mapped_orthogroup_cds.analysis.fasta",
        {
            "A": "ATGAAAACCTAA",
            "B": "ATGAAAATCTAA",
            "C": "ATGAAGACCTAA",
        },
    )
    _write(
        outdir / "hyphy" / METHOD / TRIM / "meme.json",
        json.dumps({"results": [{"site": 2, "p-value": 0.01}]}) + "\n",
    )
    _write(outdir / "asr" / METHOD / TRIM / "mlc_asr.txt", "mlc\n")
    _write(outdir / "asr" / METHOD / TRIM / "codeml_asr.ctl", "RateAncestor = 1\n")
    _write(
        outdir / "asr" / METHOD / TRIM / "rst",
        (
            "List of extant and reconstructed sequences\n"
            "node #5   ATGAAAACCTAA\n"
            "node #6   ATGAAAATCTAA\n"
            "tree with node labels for Rod Page's TreeView\n"
            "((A:0.1,B:0.1)5:0.2,C:0.1)6;\n"
        ),
    )
    return outdir


def _run_extractor(outdir: Path) -> None:
    subprocess.run(
        [
            sys.executable,
            "-m",
            "babappasnake.scripts.extract_selected_branch_ancestors",
            "--outdir",
            str(outdir),
            "--pathways",
            f"{METHOD}:{TRIM}",
            "--codeml",
            "/usr/bin/false",
            "--codeml-version-override",
            "4.10.9",
            "--gene",
            "gene1",
        ],
        check=True,
    )


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_tip_child_branch_mapping_and_sequence_recovery(tmp_path):
    outdir = _seed_common(tmp_path, "((A:0.1,B:0.1)cladeAB:0.2,C:0.1)root:0.0;", [("A", "True")])
    _run_extractor(outdir)
    rows = _read_tsv(outdir / "asr" / "branch_to_nodes.tsv")
    assert rows and rows[0]["child_node_type"] == "tip"
    assert rows[0]["status"] == "ok"


def test_internal_child_branch_mapping(tmp_path):
    outdir = _seed_common(tmp_path, "((A:0.1,B:0.1)cladeAB:0.2,C:0.1)root:0.0;", [("cladeAB", "True")])
    _run_extractor(outdir)
    rows = _read_tsv(outdir / "asr" / "branch_to_nodes.tsv")
    assert rows and rows[0]["child_node_type"] == "internal"
    assert rows[0]["status"] == "ok"


def test_multiple_selected_branches_one_pathway(tmp_path):
    outdir = _seed_common(
        tmp_path,
        "((A:0.1,B:0.1)cladeAB:0.2,C:0.1)root:0.0;",
        [("A", "True"), ("B", "True")],
    )
    _run_extractor(outdir)
    rows = _read_tsv(outdir / "asr" / "selected_branch_asr_summary.tsv")
    ok_rows = [r for r in rows if r["status"] == "ok"]
    assert len(ok_rows) >= 2


def test_no_selected_branches_skips_cleanly(tmp_path):
    outdir = _seed_common(tmp_path, "((A:0.1,B:0.1)cladeAB:0.2,C:0.1)root:0.0;", [("A", "False")])
    _run_extractor(outdir)
    rows = _read_tsv(outdir / "asr" / "selected_branch_asr_summary.tsv")
    assert rows and rows[0]["status"] == "no_selected_branches"


def test_ambiguous_mapping_failure_is_recorded(tmp_path):
    outdir = _seed_common(tmp_path, "((A-A:0.1,A_A:0.1)clX:0.2,C:0.1)root:0.0;", [("A A", "True")])
    _run_extractor(outdir)
    rows = _read_tsv(outdir / "asr" / "branch_to_nodes.tsv")
    assert rows and rows[0]["status"] == "failed"


def test_successful_extraction_writes_ancestor_descendant_fastas(tmp_path):
    outdir = _seed_common(tmp_path, "((A:0.1,B:0.1)cladeAB:0.2,C:0.1)root:0.0;", [("A", "True")])
    _run_extractor(outdir)
    anc = list(SeqIO.parse(outdir / "asr" / "ancestor_sequences_cds.fasta", "fasta"))
    desc = list(SeqIO.parse(outdir / "asr" / "descendant_sequences_cds.fasta", "fasta"))
    assert anc and desc
