from __future__ import annotations

from pathlib import Path

import pytest

from babappasnake.scripts import run_orthofinder_pipeline as of


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def test_orthofinder_stops_when_query_orthogroup_has_no_partners(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
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
    monkeypatch.setattr(
        of,
        "run_query_blast",
        lambda **kwargs: _write(
            kwargs["out_tsv"],
            "query1\tOG0001||speciesA||A1\t99\t100\t100\t100\t1e-50\t200\n",
        ),
    )
    monkeypatch.setattr(of, "make_blast_db", lambda fasta, _exe: str(Path(fasta).with_suffix("")))

    with pytest.raises(RuntimeError, match="no orthogroup partner sequences"):
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
            orthology_mode="representative",
        )


def test_orthofinder_strict_mode_stops_when_all_partners_are_multicopy(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    query = tmp_path / "query.fasta"
    proteomes = tmp_path / "proteomes"
    outdir = tmp_path / "out"
    _write(query, ">query1\nMPEPTIDE\n")
    _write(proteomes / "speciesA.fasta", ">A1\nMPEPTIDE\n>A2\nMPEPTIDE\n")

    def fake_run_cmd(cmd: list[str]) -> None:
        input_dir = Path(cmd[cmd.index("-f") + 1])
        tsv = input_dir / "OrthoFinder" / "Results_fake" / "Orthogroups" / "Orthogroups.tsv"
        _write(
            tsv,
            (
                "Orthogroup\tbabappasnake_query\tspeciesA\n"
                "OG0001\tquery1\tA1, A2\n"
            ),
        )

    monkeypatch.setattr(of, "run_cmd", fake_run_cmd)
    monkeypatch.setattr(
        of,
        "run_query_blast",
        lambda **kwargs: _write(
            kwargs["out_tsv"],
            (
                "query1\tOG0001||speciesA||A1\t99\t100\t100\t100\t1e-50\t200\n"
                "query1\tOG0001||speciesA||A2\t98\t100\t100\t100\t1e-40\t180\n"
            ),
        ),
    )
    monkeypatch.setattr(of, "make_blast_db", lambda fasta, _exe: str(Path(fasta).with_suffix("")))

    with pytest.raises(RuntimeError, match="no partner sequences retained"):
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
            orthology_mode="strict",
        )
