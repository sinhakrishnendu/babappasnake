#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import sys
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


STOP_CHARS = {'*', 'X'}
STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODON = "ATG"


def clean_protein(seq: str) -> str:
    return ''.join([aa for aa in seq if aa not in STOP_CHARS]).upper()


def trim_terminal_stop_codon(seq: str) -> str:
    seq = seq.upper()
    if len(seq) >= 3 and seq[-3:] in STOP_CODONS:
        return seq[:-3]
    return seq


def clip_lowercase_introns(raw_seq: str) -> tuple[str, bool]:
    had_lowercase = any(ch.islower() for ch in raw_seq)
    clipped = ''.join(ch for ch in raw_seq if ch.isalpha() and not ch.islower())
    return clipped.upper(), had_lowercase


def validate_full_length_cds(seq: str, rec_id: str) -> None:
    if not seq:
        raise RuntimeError(f"CDS {rec_id} is empty after lowercase intron clipping")
    if len(seq) % 3 != 0:
        raise RuntimeError(
            f"CDS {rec_id} has length {len(seq)} not divisible by 3 after lowercase intron clipping"
        )
    if not seq.startswith(START_CODON):
        raise RuntimeError(f"CDS {rec_id} does not start with uppercase start codon {START_CODON}")
    terminal = seq[-3:]
    if terminal not in STOP_CODONS:
        raise RuntimeError(
            f"CDS {rec_id} does not end with an uppercase stop codon (expected one of {sorted(STOP_CODONS)}, got {terminal})"
        )


def extract_best_orf_window(seq: str, rec_id: str) -> tuple[str, int, int]:
    starts = [i for i in range(0, len(seq) - 2) if seq[i:i + 3] == START_CODON]
    if not starts:
        raise RuntimeError(f"CDS {rec_id} has no uppercase start codon {START_CODON} after lowercase intron clipping")

    best: tuple[int, int, str] | None = None
    for start in starts:
        stop = None
        for pos in range(start + 3, len(seq) + 1, 3):
            codon = seq[pos - 3:pos]
            if codon in STOP_CODONS:
                stop = pos
                break
        if stop is None:
            continue
        candidate = seq[start:stop]
        if best is None or len(candidate) > len(best[2]):
            best = (start, stop, candidate)

    if best is None:
        raise RuntimeError(
            f"CDS {rec_id} has no in-frame uppercase stop codon downstream of an uppercase {START_CODON}"
        )

    start, stop, candidate = best
    validate_full_length_cds(candidate, rec_id)
    return candidate, start, stop


def validate_cds_for_hyphy(seq: str, rec_id: str) -> None:
    if len(seq) % 3 != 0:
        raise RuntimeError(f"Mapped CDS {rec_id} has length {len(seq)} not divisible by 3 after terminal stop trimming")
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        if codon in STOP_CODONS:
            raise RuntimeError(f"Mapped CDS {rec_id} contains internal stop codon {codon} at codon index {i // 3 + 1}")


def similarity(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    aln = pairwise2.align.globalxx(a, b, one_alignment_only=True, score_only=True)
    return aln / max(len(a), len(b))


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--proteins", required=True)
    p.add_argument("--cds", required=True)
    p.add_argument("--out-cds", required=True)
    p.add_argument("--mapping", required=True)
    p.add_argument("--out-proteins", default=None)
    p.add_argument("--min-sim", type=float, default=0.90)
    a = p.parse_args()

    protein_records = list(SeqIO.parse(a.proteins, "fasta"))
    cds_records = list(SeqIO.parse(a.cds, "fasta"))
    translated = []
    excluded = 0
    lowercase_trimmed = 0
    orf_window_clipped = 0
    for rec in cds_records:
        clipped_seq, had_lowercase = clip_lowercase_introns(str(rec.seq))
        if had_lowercase:
            lowercase_trimmed += 1
        try:
            coding_seq, start_idx, stop_idx = extract_best_orf_window(clipped_seq, rec.id)
        except RuntimeError as exc:
            excluded += 1
            print(f"[WARN] Excluding CDS {rec.id}: {exc}", file=sys.stderr)
            continue
        if start_idx != 0 or stop_idx != len(clipped_seq):
            orf_window_clipped += 1
        pep = clean_protein(str(Seq(coding_seq).translate(to_stop=False)))
        translated.append((rec, coding_seq, pep))

    if lowercase_trimmed:
        print(
            f"[INFO] Lowercase intron clipping applied to {lowercase_trimmed} CDS records.",
            file=sys.stderr,
        )
    if orf_window_clipped:
        print(
            f"[INFO] ORF window clipping applied to {orf_window_clipped} CDS records (kept uppercase ATG...STOP region).",
            file=sys.stderr,
        )
    if excluded:
        print(
            f"[INFO] Excluded {excluded} CDS records that failed start/stop or frame checks.",
            file=sys.stderr,
        )
    if not translated:
        raise RuntimeError("No valid CDS records remain after quality filtering")

    used = set()
    mapped_cds = []
    mapped_proteins = []
    rows = []
    skipped_proteins = []
    for prec in protein_records:
        pseq = clean_protein(str(prec.seq))
        best_idx = None
        best_score = -1.0
        for idx, (_, _, cpep) in enumerate(translated):
            if idx in used:
                continue
            score = similarity(pseq, cpep)
            if score > best_score:
                best_score = score
                best_idx = idx
        if best_idx is None or best_score < a.min_sim:
            skipped_proteins.append(prec.id)
            print(
                f"[WARN] Skipping protein {prec.id}: no CDS matched with translated similarity >= {a.min_sim}",
                file=sys.stderr,
            )
            continue
        used.add(best_idx)
        crec, cleaned_seq, _ = translated[best_idx]
        cds_seq = trim_terminal_stop_codon(cleaned_seq)
        validate_cds_for_hyphy(cds_seq, prec.id)
        new = SeqRecord(Seq(cds_seq), id=prec.id, description=f"mapped_from={crec.id};translated_similarity={best_score:.4f}")
        mapped_cds.append(new)
        mapped_proteins.append(SeqRecord(prec.seq, id=prec.id, description=prec.description))
        rows.append((prec.id, crec.id, f"{best_score:.4f}"))

    if skipped_proteins:
        print(
            f"[INFO] Skipped {len(skipped_proteins)} proteins with no qualifying CDS match.",
            file=sys.stderr,
        )
    if len(mapped_cds) < 3:
        raise RuntimeError(
            f"Only {len(mapped_cds)} proteins could be mapped to high-quality CDS. "
            "At least 3 are required for downstream tree-based analyses."
        )

    SeqIO.write(mapped_cds, a.out_cds, "fasta")
    if a.out_proteins:
        SeqIO.write(mapped_proteins, a.out_proteins, "fasta")
    with open(a.mapping, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["protein_header", "original_cds_header", "translated_similarity"])
        w.writerows(rows)


if __name__ == "__main__":
    main()
