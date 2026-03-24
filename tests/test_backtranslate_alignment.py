from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def test_backtranslate_alignment_is_tolerant_to_mismatch_and_short_cds(tmp_path):
    protein_aln = tmp_path / "protein.aln.fasta"
    cds_fa = tmp_path / "cds.fasta"
    out_protein = tmp_path / "out.protein.fasta"
    out_codon = tmp_path / "out.codon.fasta"
    report = tmp_path / "report.json"

    _write(
        protein_aln,
        ">seq1\nM-QK\n>seq2\nMTAK\n",
    )
    _write(
        cds_fa,
        ">seq1\nATGGCCAAATGA\n>seq2\nATGACTAAA\n",
    )

    subprocess.run(
        [
            sys.executable,
            "-m",
            "babappasnake.scripts.backtranslate_alignment",
            "--protein-aln",
            str(protein_aln),
            "--cds",
            str(cds_fa),
            "--out-protein",
            str(out_protein),
            "--out-codon",
            str(out_codon),
            "--report-json",
            str(report),
            "--min-seqs",
            "2",
        ],
        check=True,
    )

    protein_records = list(SeqIO.parse(out_protein, "fasta"))
    codon_records = list(SeqIO.parse(out_codon, "fasta"))
    assert len(protein_records) == 2
    assert len(codon_records) == 2

    codon_map = {rec.id: str(rec.seq) for rec in codon_records}
    assert len(codon_map["seq1"]) == 12
    assert len(codon_map["seq2"]) == 12
    assert "NNN" in codon_map["seq2"]

    payload = json.loads(report.read_text(encoding="utf-8"))
    seq2 = next(row for row in payload["records"] if row["protein_id"] == "seq2")
    assert int(seq2["padded_short_codons"]) >= 1
