#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


STOP_CODONS = {"TAA", "TAG", "TGA"}


def strip_terminal_stop_in_aligned_seq(aligned_seq: str, rec_id: str) -> tuple[str, bool]:
    chars = list(aligned_seq.upper())
    nongap_idx = [i for i, c in enumerate(chars) if c != "-"]
    ungapped = "".join(chars[i] for i in nongap_idx)

    if len(ungapped) % 3 != 0:
        raise RuntimeError(
            f"Sequence {rec_id} has ungapped length {len(ungapped)} not divisible by 3; cannot process codon alignment safely"
        )

    if len(ungapped) >= 3 and ungapped[-3:] in STOP_CODONS:
        for idx in nongap_idx[-3:]:
            chars[idx] = "-"
        return "".join(chars), True

    return "".join(chars), False


def drop_all_gap_codon_columns(seqs: list[str]) -> tuple[list[str], int]:
    if not seqs:
        return seqs, 0

    aln_len = len(seqs[0])
    if any(len(s) != aln_len for s in seqs):
        raise RuntimeError("Input alignment sequences do not all have the same length")
    if aln_len % 3 != 0:
        raise RuntimeError(f"Aligned length {aln_len} is not divisible by 3")

    kept_codon_starts = []
    removed = 0
    for i in range(0, aln_len, 3):
        if all(s[i : i + 3] == "---" for s in seqs):
            removed += 1
            continue
        kept_codon_starts.append(i)

    trimmed = ["".join(s[i : i + 3] for i in kept_codon_starts) for s in seqs]
    return trimmed, removed


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    a = p.parse_args()

    records = list(SeqIO.parse(a.input, "fasta"))
    if not records:
        raise RuntimeError(f"No sequences found in {a.input}")

    updated_records = []
    stripped_count = 0
    seqs = []
    for rec in records:
        new_seq, stripped = strip_terminal_stop_in_aligned_seq(str(rec.seq), rec.id)
        stripped_count += int(stripped)
        seqs.append(new_seq)
        updated_records.append(
            SeqRecord(Seq(new_seq), id=rec.id, description=rec.description)
        )

    trimmed_seqs, removed_codon_cols = drop_all_gap_codon_columns(seqs)
    for rec, seq in zip(updated_records, trimmed_seqs):
        rec.seq = Seq(seq)

    out_path = Path(a.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(updated_records, str(out_path), "fasta")

    print(
        f"Wrote {len(updated_records)} sequences to {out_path} "
        f"(terminal stops stripped: {stripped_count}; all-gap codon columns removed: {removed_codon_cols})"
    )


if __name__ == "__main__":
    main()
