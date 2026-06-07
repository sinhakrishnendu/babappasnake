from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from babappasnake.scripts.prepare_external_orthogroup import prepare_external_orthogroup


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def test_prepare_external_orthogroup_writes_standard_outputs(tmp_path: Path):
    proteins = tmp_path / "external.fasta"
    query = tmp_path / "query.fasta"
    outdir = tmp_path / "out"
    _write(proteins, ">q1\nMPEPTIDE\n>s1\nMPEPTIDE\n>s2\nMPEPTIDE\n")
    _write(query, ">q1\nMPEPTIDE\n")

    prepare_external_orthogroup(proteins, outdir, query_fasta=str(query))

    assert (outdir / "orthogroup_proteins.fasta").exists()
    assert (outdir / "orthogroup_headers.txt").read_text(encoding="utf-8").splitlines() == ["q1", "s1", "s2"]
    rows = list(csv.DictReader((outdir / "orthogroup_summary.tsv").open(encoding="utf-8"), delimiter="\t"))
    assert [row["selected_members"] for row in rows] == ["s1", "s2"]
    assert {row["orthology_mode"] for row in rows} == {"external"}
    metadata = json.loads((outdir / "orthogroup_metadata.json").read_text(encoding="utf-8"))
    assert metadata["source"] == "external"
    assert metadata["retained_partner_count"] == 2


def test_prepare_external_orthogroup_rejects_duplicate_ids(tmp_path: Path):
    proteins = tmp_path / "external.fasta"
    _write(proteins, ">q1\nMPEPTIDE\n>q1\nMPEPTIDE\n")

    with pytest.raises(RuntimeError, match="duplicate sequence IDs"):
        prepare_external_orthogroup(proteins, tmp_path / "out")
