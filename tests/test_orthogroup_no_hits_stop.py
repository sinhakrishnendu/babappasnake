from __future__ import annotations

from pathlib import Path

import pytest

from babappasnake.scripts import run_orthofinder_pipeline as of
from babappasnake.scripts import run_rbh_pipeline as rbh


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def test_rbh_fallback_hard_stops_when_both_methods_have_zero_one_to_one(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    query = tmp_path / "query.fasta"
    proteomes = tmp_path / "proteomes"
    outdir = tmp_path / "out"
    _write(query, ">query1\nMPEPTIDE\n")
    _write(proteomes / "speciesA.fasta", ">A1\nMPEPTIDE\n")

    monkeypatch.setattr(rbh, "make_blast_db", lambda fasta, exe: str(Path(fasta).with_suffix("")))
    monkeypatch.setattr(rbh, "run_blastp", lambda *args, **kwargs: None)
    monkeypatch.setattr(rbh, "get_best_hits", lambda *args, **kwargs: {})

    def fake_fallback(*args, **kwargs):
        fallback_dir = kwargs["outdir"]
        _write(
            fallback_dir / "rbh_summary.tsv",
            "species\tquery\tortholog\nspeciesA\tquery1\tNA\n",
        )
        _write(fallback_dir / "orthogroup_headers.txt", "query1\n")
        _write(fallback_dir / "orthogroup_proteins.fasta", ">query1\nMPEPTIDE\n")

    monkeypatch.setattr(rbh, "run_orthofinder_fallback", fake_fallback)

    with pytest.raises(RuntimeError, match="after RBH and OrthoFinder fallback"):
        rbh.choose_single_copy_rbh(
            query_fasta=query,
            proteomes_dir=proteomes,
            outdir=outdir,
            coverage=0.7,
            evalue=1e-5,
            threads=1,
            blastp_exe="blastp",
            makeblastdb_exe="makeblastdb",
            orthofinder_exe="orthofinder",
        )


def test_orthofinder_stops_when_query_orthogroup_has_no_partners(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    query = tmp_path / "query.fasta"
    proteomes = tmp_path / "proteomes"
    outdir = tmp_path / "out"
    _write(query, ">query1\nMPEPTIDE\n")
    _write(proteomes / "speciesA.fasta", ">A1\nMPEPTIDE\n")

    def fake_run_cmd(cmd: list[str]) -> None:
        input_dir = Path(cmd[cmd.index("-f") + 1])
        tsv = input_dir / "OrthoFinder" / "Results_fake" / "Orthogroups" / "Orthogroups.tsv"
        _write(
            tsv,
            (
                "Orthogroup\tbabappasnake_query\tspeciesA\n"
                "OG0001\tquery1\t\n"
            ),
        )

    monkeypatch.setattr(of, "run_cmd", fake_run_cmd)
    monkeypatch.setattr(of, "run_query_blast", lambda **kwargs: _write(
        kwargs["out_tsv"],
        "query1\tOG0001||speciesA||A1\t99\t100\t100\t100\t1e-50\t200\n",
    ))
    monkeypatch.setattr(of, "make_blast_db", lambda fasta, _exe: str(Path(fasta).with_suffix("")))

    with pytest.raises(RuntimeError, match="no orthogroup partner sequences for BLAST-based query mapping"):
        of.choose_orthogroup_with_orthofinder(
            query_fasta=query,
            proteomes_dir=proteomes,
            outdir=outdir,
            threads=1,
            coverage=0.7,
            evalue=1e-5,
            blastp_exe="blastp",
            makeblastdb_exe="makeblastdb",
            orthofinder_exe="orthofinder",
        )


def test_rbh_fallback_selects_orthofinder_when_it_has_more_one_to_one(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    query = tmp_path / "query.fasta"
    proteomes = tmp_path / "proteomes"
    outdir = tmp_path / "out"
    _write(query, ">query1\nMPEPTIDE\n")
    _write(proteomes / "speciesA.fasta", ">A1\nMPEPTIDE\n")
    _write(proteomes / "speciesB.fasta", ">B1\nMPEPTIDE\n")

    monkeypatch.setattr(rbh, "make_blast_db", lambda fasta, exe: str(Path(fasta).with_suffix("")))
    monkeypatch.setattr(rbh, "run_blastp", lambda *args, **kwargs: None)

    def fake_best_hits(blast_file: Path, _min_cov: float):
        name = blast_file.name
        if "speciesA" in name:
            if name.endswith(".fwd.tsv"):
                return {"query1": "A1"}
            return {"A1": "query1"}
        return {}

    monkeypatch.setattr(rbh, "get_best_hits", fake_best_hits)

    def fake_fallback(*args, **kwargs):
        fallback_dir = kwargs["outdir"]
        _write(
            fallback_dir / "rbh_summary.tsv",
            (
                "species\tquery\tortholog\n"
                "speciesA\tquery1\tA1\n"
                "speciesB\tquery1\tB1\n"
            ),
        )
        _write(fallback_dir / "orthogroup_headers.txt", "query1\nA1\nB1\n")
        _write(
            fallback_dir / "orthogroup_proteins.fasta",
            ">query1\nMPEPTIDE\n>A1\nMPEPTIDE\n>B1\nMPEPTIDE\n",
        )

    monkeypatch.setattr(rbh, "run_orthofinder_fallback", fake_fallback)

    rbh.choose_single_copy_rbh(
        query_fasta=query,
        proteomes_dir=proteomes,
        outdir=outdir,
        coverage=0.7,
        evalue=1e-5,
        threads=1,
        blastp_exe="blastp",
        makeblastdb_exe="makeblastdb",
        orthofinder_exe="orthofinder",
    )

    summary = (outdir / "rbh_summary.tsv").read_text(encoding="utf-8")
    assert "speciesA\tquery1\tA1" in summary
    assert "speciesB\tquery1\tB1" in summary
