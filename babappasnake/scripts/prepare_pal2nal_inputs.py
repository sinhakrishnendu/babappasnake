#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO
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

    for rec in protein_records:
        prot_id = rec.id.strip()
        if not prot_id:
            raise RuntimeError("Encountered empty protein ID in alignment")
        if prot_id in seen_ids:
            raise RuntimeError(f"Duplicate protein ID in alignment: {prot_id}")
        seen_ids.add(prot_id)

        cds_rec = resolve_cds_for_protein(rec, cds_index)
        normalized_protein.append(SeqRecord(rec.seq, id=prot_id, description=""))
        ordered_cds.append(SeqRecord(cds_rec.seq, id=prot_id, description=""))

    out_protein = Path(a.out_protein)
    out_cds = Path(a.out_cds)
    out_protein.parent.mkdir(parents=True, exist_ok=True)
    out_cds.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(normalized_protein, out_protein, "fasta")
    SeqIO.write(ordered_cds, out_cds, "fasta")


if __name__ == "__main__":
    main()

