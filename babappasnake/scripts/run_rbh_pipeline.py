#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path
from typing import Dict

from Bio import SeqIO


FASTA_EXTENSIONS = {".fa", ".faa", ".fasta"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--query", required=True)
    p.add_argument("--proteomes", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--coverage", type=float, default=0.70)
    p.add_argument("--evalue", type=float, default=1e-5)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--blastp", default="blastp")
    p.add_argument("--makeblastdb", default="makeblastdb")
    return p.parse_args()


def run_cmd(cmd: list[str]) -> None:
    result = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed:\n{' '.join(cmd)}\n\nSTDERR:\n{result.stderr}")


def _warn_skip(path: Path, reason: str) -> None:
    print(f"[WARN] Skipping proteome file '{path}': {reason}", file=sys.stderr)


def _fasta_sanity_issue(path: Path) -> str | None:
    try:
        if path.stat().st_size == 0:
            return "empty file"
    except OSError as exc:
        return f"unreadable file ({exc})"

    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    return None
                return "first non-empty line does not start with '>'"
    except OSError as exc:
        return f"unreadable file ({exc})"

    return "no non-empty lines"


def looks_like_fasta(path: Path) -> bool:
    return _fasta_sanity_issue(path) is None


def discover_proteome_fastas(proteomes_dir: Path) -> list[Path]:
    if not proteomes_dir.is_dir():
        raise FileNotFoundError(f"Proteomes directory not found: {proteomes_dir}")

    discovered: list[Path] = []
    for path in sorted(proteomes_dir.iterdir()):
        name = path.name
        if not path.is_file():
            continue
        if name.startswith(".") or name.startswith("._") or name == ".DS_Store":
            _warn_skip(path, "hidden/macOS metadata file")
            continue
        if path.suffix.lower() not in FASTA_EXTENSIONS:
            _warn_skip(path, f"unsupported extension '{path.suffix}'")
            continue
        issue = _fasta_sanity_issue(path)
        if issue is not None:
            _warn_skip(path, issue)
            continue
        discovered.append(path)

    if not discovered:
        raise FileNotFoundError(
            f"No valid proteome FASTA files found in {proteomes_dir}. "
            "Expected .fa/.faa/.fasta files with FASTA headers and non-empty content."
        )
    return discovered


def make_blast_db(fasta: Path, makeblastdb_exe: str) -> str:
    prefix = str(fasta.with_suffix(""))
    if not Path(prefix + ".pin").exists():
        run_cmd([makeblastdb_exe, "-in", str(fasta), "-dbtype", "prot", "-out", prefix])
    return prefix


def run_blastp(query: Path, db_prefix: str, output: Path, evalue: float, threads: int, blastp_exe: str) -> None:
    cmd = [
        blastp_exe, "-query", str(query), "-db", db_prefix,
        "-out", str(output), "-evalue", str(evalue),
        "-outfmt", "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore",
        "-num_threads", str(threads), "-seg", "yes", "-soft_masking", "true", "-use_sw_tback",
    ]
    run_cmd(cmd)


def get_best_hits(blast_file: Path, min_cov: float) -> Dict[str, str]:
    best = {}
    with open(blast_file, encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            q, s, _pident, length, qlen, slen, *_rest, _evalue, bitscore = parts
            length = int(length)
            qlen = int(qlen)
            slen = int(slen)
            bitscore = float(bitscore)
            qcov = length / max(1, qlen)
            scov = length / max(1, slen)
            if qcov >= min_cov and scov >= min_cov:
                if q not in best or bitscore > best[q][1]:
                    best[q] = (s, bitscore)
    return {q: s for q, (s, _b) in best.items()}


def reciprocal_best_hits(fwd: Dict[str, str], rev: Dict[str, str]) -> Dict[str, str]:
    return {q: s for q, s in fwd.items() if s in rev and rev[s] == q}


def choose_single_copy_rbh(query_fasta: Path, proteomes_dir: Path, outdir: Path, coverage: float, evalue: float, threads: int, blastp_exe: str, makeblastdb_exe: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    qdb = make_blast_db(query_fasta, makeblastdb_exe)
    species_fastas = discover_proteome_fastas(proteomes_dir)

    query_records = list(SeqIO.parse(str(query_fasta), "fasta"))
    if len(query_records) != 1:
        raise ValueError("This simple version expects exactly one query protein in query.fasta")
    query_id = query_records[0].id

    summary_rows = []
    ortholog_records = [query_records[0]]
    headers = [query_id]

    for proteome in species_fastas:
        species = proteome.stem
        sp_db = make_blast_db(proteome, makeblastdb_exe)
        fwd = outdir / f"{species}.fwd.tsv"
        rev = outdir / f"{species}.rev.tsv"
        run_blastp(query_fasta, sp_db, fwd, evalue, threads, blastp_exe)
        run_blastp(proteome, qdb, rev, evalue, threads, blastp_exe)
        fwd_hits = get_best_hits(fwd, coverage)
        rev_hits = get_best_hits(rev, coverage)
        rbh = reciprocal_best_hits(fwd_hits, rev_hits)
        target = rbh.get(query_id)
        summary_rows.append((species, query_id, target if target else "NA"))
        if target:
            wanted = None
            for rec in SeqIO.parse(str(proteome), "fasta"):
                if rec.id == target:
                    wanted = rec
                    break
            if wanted is None:
                raise RuntimeError(f"RBH header {target} not found in {proteome}")
            ortholog_records.append(wanted)
            headers.append(wanted.id)

    SeqIO.write(ortholog_records, outdir / "orthogroup_proteins.fasta", "fasta")
    with open(outdir / "orthogroup_headers.txt", "w", encoding="utf-8") as fh:
        for h in headers:
            fh.write(h + "\n")
    with open(outdir / "rbh_summary.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["species", "query", "ortholog"])
        w.writerows(summary_rows)


if __name__ == "__main__":
    a = parse_args()
    choose_single_copy_rbh(Path(a.query), Path(a.proteomes), Path(a.outdir), a.coverage, a.evalue, a.threads, a.blastp, a.makeblastdb)
