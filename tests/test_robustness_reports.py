from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path


METHODS = ("babappalign", "mafft", "prank")
TRIM_STATES = ("raw", "clipkit")


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _write_fasta(path: Path, records: dict[str, str]) -> None:
    lines = []
    for name, seq in records.items():
        lines.append(f">{name}")
        lines.append(seq)
    _write(path, "\n".join(lines) + "\n")


def _seed_pathway(outdir: Path, method: str, trim_state: str) -> None:
    _write_fasta(
        outdir / "alignments" / method / trim_state / "orthogroup_proteins.analysis.fasta",
        {"seqA": "MKTAA", "seqB": "MKTEA"},
    )
    _write_fasta(
        outdir / "alignments" / method / trim_state / "mapped_orthogroup_cds.analysis.fasta",
        {"seqA": "ATGAAAACCGCTTAA", "seqB": "ATGAAAACCGAATAA"},
    )
    _write(outdir / "tree" / method / trim_state / "orthogroup.rooted.treefile", "(seqA,seqB);\n")

    _write(
        outdir / "hyphy" / method / trim_state / "significant_foregrounds.tsv",
        "foreground_branch\tabsrel_p\nbranchA\t0.01\n",
    )
    _write(
        outdir / "hyphy" / method / trim_state / "foreground_threshold.json",
        json.dumps({"selected_pcut": 0.05}) + "\n",
    )
    _write(
        outdir / "hyphy" / method / trim_state / "absrel.json",
        json.dumps(
            {
                "branch attributes": {
                    "0": {
                        "branchA": {"is leaf": True},
                        "internal1": {"is leaf": False},
                    }
                }
            }
        )
        + "\n",
    )
    _write(
        outdir / "hyphy" / method / trim_state / "meme.json",
        json.dumps({"results": [{"site": "10", "p-value": 0.01}]}) + "\n",
    )

    _write(
        outdir / "branchsite" / method / trim_state / "branchsite_results.tsv",
        (
            "foreground_branch\tLRT\tp_value\tq_value\tsignificant_BH_0.05\n"
            "branchA\t12.3\t0.001\t0.01\tTrue\n"
        ),
    )
    _write(
        outdir / "hyphy" / method / trim_state / "hyphy_done.json",
        json.dumps({"absrel_status": "ok", "meme_status": "ok"}) + "\n",
    )
    _write(outdir / "asr" / method / trim_state / "asr_done.json", '{"done": true}\n')


def test_generate_robustness_reports_for_six_method_trim_pathways(tmp_path):
    outdir = tmp_path / "run"
    for method in METHODS:
        for trim_state in TRIM_STATES:
            _seed_pathway(outdir, method, trim_state)

    summary_dir = outdir / "summary"
    matrix_out = summary_dir / "robustness_matrix.tsv"
    consensus_out = summary_dir / "robustness_consensus.tsv"
    narrative_out = summary_dir / "robustness_narrative.txt"
    comparative_out = summary_dir / "comparative_reproducibility_summary.txt"
    latex_out = summary_dir / "robustness_publication_table.tex"
    pathways = ",".join(f"{m}:{t}" for m in METHODS for t in TRIM_STATES)

    subprocess.run(
        [
            sys.executable,
            "-m",
            "babappasnake.scripts.generate_robustness_reports",
            "--outdir",
            str(outdir),
            "--pathways",
            pathways,
            "--meme-p",
            "0.05",
            "--matrix-out",
            str(matrix_out),
            "--consensus-out",
            str(consensus_out),
            "--narrative-out",
            str(narrative_out),
            "--comparative-out",
            str(comparative_out),
            "--latex-out",
            str(latex_out),
        ],
        check=True,
    )

    assert matrix_out.exists()
    assert consensus_out.exists()
    assert narrative_out.exists()
    assert comparative_out.exists()
    assert latex_out.exists()

    with open(matrix_out, encoding="utf-8") as fh:
        matrix_rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(matrix_rows) == 6
    assert len({(row["method"], row["trim_state"]) for row in matrix_rows}) == 6
    assert len({row["pathway"] for row in matrix_rows}) == 6
    assert all(row["status"] == "ok" for row in matrix_rows)

    with open(consensus_out, encoding="utf-8") as fh:
        consensus_rows = list(csv.DictReader(fh, delimiter="\t"))

    branch_row = next(
        row
        for row in consensus_rows
        if row["signal_source"] == "aBSREL_branch" and row["signal_id"] == "branchA"
    )
    assert branch_row["replication_count"] == "6"
    assert branch_row["interpretation"] == "highly_robust"

    comparative_text = comparative_out.read_text(encoding="utf-8")
    assert "Pathways evaluated (6)" in comparative_text

    latex_text = latex_out.read_text(encoding="utf-8")
    assert r"aBSREL\_branch" in latex_text
    assert r"aBSREL\\_branch" not in latex_text
