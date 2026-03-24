#!/usr/bin/env python3
from __future__ import annotations

import argparse
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
    keys = []
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


def resolve_cds_for_protein(rec: SeqRecord, cds_index: dict[str, SeqRecord]) -> SeqRecord:
    for key in collect_record_keys(rec):
        hit = cds_index.get(key)
        if hit is not None:
            return hit
    raise RuntimeError(
        f"No CDS match found for aligned protein header '{rec.description or rec.id}'. "
        "Run bl2seq (-p tblastn) or GeneWise to inspect inconsistencies."
    )


def translated_cds_protein(cds_rec: SeqRecord) -> str:
    translated = str(Seq(str(cds_rec.seq)).translate(to_stop=False)).upper()
    if translated.endswith("*"):
        translated = translated[:-1]
    return translated


def reconcile_alignment_to_translation(aligned_protein: str, translated: str) -> str:
    nongap = [ch for ch in aligned_protein if ch != "-"]
    if len(nongap) != len(translated):
        raise RuntimeError(
            f"Alignment residue count ({len(nongap)}) does not match CDS-translated length ({len(translated)})"
        )
    out: list[str] = []
    idx = 0
    for ch in aligned_protein:
        if ch == "-":
            out.append("-")
        else:
            out.append(translated[idx])
            idx += 1
    return "".join(out)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--protein-aln", required=True)
    p.add_argument("--cds", required=True)
    p.add_argument("--out-protein", required=True)
    p.add_argument("--out-cds", required=True)
    a = p.parse_args()

    protein_records = list(SeqIO.parse(a.protein_aln, "fasta"))
    cds_records = list(SeqIO.parse(a.cds, "fasta"))
    if not protein_records:
        raise RuntimeError("Protein alignment is empty")
    if not cds_records:
        raise RuntimeError("Mapped CDS FASTA is empty")

    cds_index = build_cds_index(cds_records)
    normalized_protein: list[SeqRecord] = []
    ordered_cds: list[SeqRecord] = []
    seen_ids: set[str] = set()
    dropped: list[str] = []

    for rec in protein_records:
        prot_id = rec.id.strip()
        if not prot_id:
            raise RuntimeError("Encountered empty protein ID in alignment")
        if prot_id in seen_ids:
            raise RuntimeError(f"Duplicate protein ID in alignment: {prot_id}")
        seen_ids.add(prot_id)

        cds_rec = resolve_cds_for_protein(rec, cds_index)
        translated = translated_cds_protein(cds_rec)
        try:
            reconciled = reconcile_alignment_to_translation(str(rec.seq).upper(), translated)
        except RuntimeError as exc:
            dropped.append(f"{prot_id}: {exc}")
            continue
        normalized_protein.append(SeqRecord(Seq(reconciled), id=prot_id, description=""))
        ordered_cds.append(SeqRecord(cds_rec.seq, id=prot_id, description=""))

    if dropped:
        print(
            f"[WARN] Dropped {len(dropped)} sequence(s) due protein/CDS length inconsistency before pal2nal.",
            file=sys.stderr,
        )
        for msg in dropped:
            print(f"[WARN] {msg}", file=sys.stderr)

    if len(normalized_protein) < 3:
        raise RuntimeError(
            f"Only {len(normalized_protein)} sequence(s) remained after pal2nal consistency filtering; at least 3 are required."
        )

    out_protein = Path(a.out_protein)
    out_cds = Path(a.out_cds)
    out_protein.parent.mkdir(parents=True, exist_ok=True)
    out_cds.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(normalized_protein, out_protein, "fasta")
    SeqIO.write(ordered_cds, out_cds, "fasta")


if __name__ == "__main__":
    main()
