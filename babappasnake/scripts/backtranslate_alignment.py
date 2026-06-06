#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def header_keys(text: str) -> list[str]:
    raw = str(text).strip()
    if not raw:
        return []
    first = raw.split()[0]
    keys: list[str] = []
    for key in (raw, first):
        if key and key not in keys:
            keys.append(key)
    return keys


def collect_record_keys(rec: SeqRecord) -> list[str]:
    keys: list[str] = []
    for source in (rec.id, rec.description):
        for key in header_keys(source):
            if key not in keys:
                keys.append(key)
    return keys


def build_cds_index(records: list[SeqRecord]) -> dict[str, SeqRecord]:
    index: dict[str, SeqRecord] = {}
    for rec in records:
        for key in collect_record_keys(rec):
            if key in index and index[key].id != rec.id:
                raise RuntimeError(
                    f"Ambiguous CDS header key '{key}' maps to both {index[key].id} and {rec.id}"
                )
            index[key] = rec
    return index


def clean_cds(seq: str) -> str:
    letters = [ch for ch in str(seq).upper() if ch.isalpha()]
    return "".join(letters)


def codon_to_aa(codon: str) -> str:
    if len(codon) != 3:
        return "X"
    if any(base not in {"A", "C", "G", "T"} for base in codon):
        return "X"
    aa = str(Seq(codon).translate(to_stop=False)).upper()
    return aa[0] if aa else "X"


def reconcile_record(
    protein_aln_seq: str,
    cds_seq: str,
) -> tuple[str, str, dict[str, int]]:
    aligned = str(protein_aln_seq).upper()
    cleaned_cds = clean_cds(cds_seq)
    full_codons = len(cleaned_cds) // 3
    codons = [cleaned_cds[i * 3:(i + 1) * 3] for i in range(full_codons)]
    trailing_nt = len(cleaned_cds) % 3

    codon_idx = 0
    pad_short = 0
    aa_mismatch = 0
    internal_stop = 0
    nongap = 0

    protein_out: list[str] = []
    codon_out: list[str] = []
    for aa_char in aligned:
        if aa_char == "-":
            protein_out.append("-")
            codon_out.append("---")
            continue

        nongap += 1
        if codon_idx < len(codons):
            codon = codons[codon_idx]
            codon_idx += 1
        else:
            codon = "NNN"
            pad_short += 1

        aa_from_codon = codon_to_aa(codon)
        if aa_from_codon == "*" and codon_idx < len(codons):
            internal_stop += 1

        expected = aa_char.upper()
        if expected not in {"X", "?", "*", "B", "Z", "J", "U", "O"} and aa_from_codon != expected:
            aa_mismatch += 1

        protein_out.append(aa_from_codon if aa_from_codon != "*" else "X")
        codon_out.append(codon)

    leftover_codons = max(0, len(codons) - codon_idx)
    stats = {
        "nongap_residues": nongap,
        "input_codon_count": len(codons),
        "padded_short_codons": pad_short,
        "leftover_codons": leftover_codons,
        "aa_mismatches": aa_mismatch,
        "internal_stops": internal_stop,
        "trailing_nt_trimmed": trailing_nt,
    }
    return "".join(protein_out), "".join(codon_out), stats


def main() -> None:
    p = argparse.ArgumentParser(
        description=(
            "Robust pal2nal-like back-translation: projects CDS onto an aligned protein MSA "
            "while preserving sequence order and IDs and tolerating mismatches."
        )
    )
    p.add_argument("--protein-aln", required=True)
    p.add_argument("--cds", required=True)
    p.add_argument("--out-protein", required=True)
    p.add_argument("--out-codon", required=True)
    p.add_argument("--report-json", required=False, default="")
    p.add_argument("--min-seqs", type=int, default=3)
    a = p.parse_args()

    protein_records = list(SeqIO.parse(a.protein_aln, "fasta"))
    cds_records = list(SeqIO.parse(a.cds, "fasta"))
    if not protein_records:
        raise RuntimeError("Protein alignment is empty.")
    if not cds_records:
        raise RuntimeError("CDS FASTA is empty.")

    cds_index = build_cds_index(cds_records)
    seen_ids: set[str] = set()
    kept_protein: list[SeqRecord] = []
    kept_codon: list[SeqRecord] = []

    dropped_no_cds = 0
    report_rows: list[dict[str, object]] = []

    for rec in protein_records:
        prot_id = rec.id.strip()
        if not prot_id:
            raise RuntimeError("Encountered empty protein ID in protein alignment.")
        if prot_id in seen_ids:
            raise RuntimeError(f"Duplicate protein ID in protein alignment: {prot_id}")
        seen_ids.add(prot_id)

        cds_hit = None
        for key in collect_record_keys(rec):
            cds_hit = cds_index.get(key)
            if cds_hit is not None:
                break

        if cds_hit is None:
            dropped_no_cds += 1
            print(
                f"[WARN] Skipping protein '{prot_id}' because no CDS header match was found.",
                file=sys.stderr,
            )
            continue

        reconciled_protein, codon_aln, stats = reconcile_record(str(rec.seq), str(cds_hit.seq))
        kept_protein.append(SeqRecord(Seq(reconciled_protein), id=prot_id, description=""))
        kept_codon.append(SeqRecord(Seq(codon_aln), id=prot_id, description=""))
        report_rows.append({"protein_id": prot_id, **stats})

    if dropped_no_cds:
        print(f"[WARN] Dropped {dropped_no_cds} sequence(s) due missing CDS matches.", file=sys.stderr)

    if len(kept_protein) < int(a.min_seqs):
        raise RuntimeError(
            f"Only {len(kept_protein)} sequence(s) remained after robust back-translation; "
            f"at least {a.min_seqs} are required."
        )

    out_protein = Path(a.out_protein)
    out_codon = Path(a.out_codon)
    out_protein.parent.mkdir(parents=True, exist_ok=True)
    out_codon.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(kept_protein, out_protein, "fasta")
    SeqIO.write(kept_codon, out_codon, "fasta")

    if a.report_json:
        payload = {
            "n_input_proteins": len(protein_records),
            "n_input_cds": len(cds_records),
            "n_kept": len(kept_protein),
            "n_dropped_no_cds": dropped_no_cds,
            "records": report_rows,
        }
        report_path = Path(a.report_json)
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
