from __future__ import annotations

from pathlib import Path

import pytest

from babappasnake.scripts.run_rbh_pipeline import discover_proteome_fastas


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def test_valid_fasta_accepted(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    _write(tmp_path / "speciesA.fasta", ">a\nMPEPTIDE\n")
    files = discover_proteome_fastas(tmp_path)
    assert [p.name for p in files] == ["speciesA.fasta"]
    assert capsys.readouterr().err == ""


def test_appledouble_sidecar_rejected(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    _write(tmp_path / "speciesA.fasta", ">a\nMPEPTIDE\n")
    _write(tmp_path / "._speciesA.fasta", "not fasta\n")
    files = discover_proteome_fastas(tmp_path)
    assert [p.name for p in files] == ["speciesA.fasta"]
    assert "._speciesA.fasta" in capsys.readouterr().err


def test_ds_store_rejected(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    _write(tmp_path / "speciesA.faa", ">a\nMPEPTIDE\n")
    _write(tmp_path / ".DS_Store", "junk\n")
    files = discover_proteome_fastas(tmp_path)
    assert [p.name for p in files] == ["speciesA.faa"]
    assert ".DS_Store" in capsys.readouterr().err


def test_hidden_file_rejected(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    _write(tmp_path / "speciesA.fa", ">a\nMPEPTIDE\n")
    _write(tmp_path / ".hidden.fa", ">x\nMXXX\n")
    files = discover_proteome_fastas(tmp_path)
    assert [p.name for p in files] == ["speciesA.fa"]
    assert ".hidden.fa" in capsys.readouterr().err


def test_malformed_fasta_rejected(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    _write(tmp_path / "speciesA.fasta", ">a\nMPEPTIDE\n")
    _write(tmp_path / "bad.fasta", "MPEPTIDE_WITHOUT_HEADER\n")
    files = discover_proteome_fastas(tmp_path)
    assert [p.name for p in files] == ["speciesA.fasta"]
    err = capsys.readouterr().err
    assert "bad.fasta" in err
    assert "does not start with '>'" in err


def test_raises_when_no_valid_fasta_found(tmp_path: Path) -> None:
    _write(tmp_path / "._species.fasta", "junk\n")
    _write(tmp_path / ".DS_Store", "junk\n")
    _write(tmp_path / "bad.fasta", "MPEPTIDE\n")
    with pytest.raises(FileNotFoundError, match="No valid proteome FASTA files found"):
        discover_proteome_fastas(tmp_path)
