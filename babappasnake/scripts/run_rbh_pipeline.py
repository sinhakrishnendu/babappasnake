#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path
from typing import Dict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


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
    p.add_argument("--orthofinder", default="orthofinder")
    p.add_argument("--mode", choices=["rbh", "rbh_fallback"], default="rbh")
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


def write_orthogroup_outputs(
    outdir: Path,
    ortholog_records: list[SeqRecord],
    headers: list[str],
    summary_rows: list[tuple[str, str, str]],
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    SeqIO.write(ortholog_records, outdir / "orthogroup_proteins.fasta", "fasta")
    with open(outdir / "orthogroup_headers.txt", "w", encoding="utf-8") as fh:
        for header in headers:
            fh.write(header + "\n")
    with open(outdir / "rbh_summary.tsv", "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["species", "query", "ortholog"])
        writer.writerows(summary_rows)


def count_one_to_one_orthologs(summary_rows: list[tuple[str, str, str]]) -> int:
    count = 0
    for _species, _query, ortholog in summary_rows:
        value = str(ortholog or "").strip()
        if not value or value == "NA":
            continue
        members = [token.strip() for token in value.split(",") if token.strip()]
        if len(members) == 1:
            count += 1
    return count


def read_summary_rows(summary_tsv: Path) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    with open(summary_tsv, encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(
                (
                    str(row.get("species", "")).strip(),
                    str(row.get("query", "")).strip(),
                    str(row.get("ortholog", "")).strip(),
                )
            )
    return rows


def run_orthofinder_fallback(
    query_fasta: Path,
    proteomes_dir: Path,
    outdir: Path,
    coverage: float,
    evalue: float,
    threads: int,
    blastp_exe: str,
    makeblastdb_exe: str,
    orthofinder_exe: str,
) -> None:
    cmd = [
        sys.executable,
        "-m",
        "babappasnake.scripts.run_orthofinder_pipeline",
        "--query",
        str(query_fasta),
        "--proteomes",
        str(proteomes_dir),
        "--outdir",
        str(outdir),
        "--coverage",
        str(float(coverage)),
        "--evalue",
        str(float(evalue)),
        "--threads",
        str(max(1, int(threads))),
        "--blastp",
        str(blastp_exe),
        "--makeblastdb",
        str(makeblastdb_exe),
        "--orthofinder",
        str(orthofinder_exe),
    ]
    run_cmd(cmd)


def run_rbh_only(
    query_fasta: Path,
    proteomes_dir: Path,
    outdir: Path,
    coverage: float,
    evalue: float,
    threads: int,
    blastp_exe: str,
    makeblastdb_exe: str,
) -> tuple[str, list[SeqRecord], list[str], list[tuple[str, str, str]], int]:
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
    return query_id, ortholog_records, headers, summary_rows, len(species_fastas)


def choose_single_copy_rbh(
    query_fasta: Path,
    proteomes_dir: Path,
    outdir: Path,
    coverage: float,
    evalue: float,
    threads: int,
    blastp_exe: str,
    makeblastdb_exe: str,
    orthofinder_exe: str,
    mode: str = "rbh",
) -> None:
    selected_mode = str(mode).strip().lower()
    if selected_mode not in {"rbh", "rbh_fallback"}:
        raise ValueError(f"Unsupported RBH selection mode: {mode}")

    outdir.mkdir(parents=True, exist_ok=True)
    query_id, rbh_records, rbh_headers, rbh_rows, total_species = run_rbh_only(
        query_fasta=query_fasta,
        proteomes_dir=proteomes_dir,
        outdir=outdir,
        coverage=coverage,
        evalue=evalue,
        threads=threads,
        blastp_exe=blastp_exe,
        makeblastdb_exe=makeblastdb_exe,
    )
    write_orthogroup_outputs(outdir=outdir, ortholog_records=rbh_records, headers=rbh_headers, summary_rows=rbh_rows)

    rbh_one_to_one = count_one_to_one_orthologs(rbh_rows)
    total = max(1, int(total_species))
    rbh_ratio = float(rbh_one_to_one) / float(total)
    print(
        f"[INFO] RBH 1:1 ortholog recovery for query '{query_id}': "
        f"{rbh_one_to_one}/{total} ({rbh_ratio:.1%})."
    )

    if selected_mode == "rbh":
        if rbh_one_to_one <= 0:
            raise RuntimeError(
                f"No orthogroup could be defined for query '{query_id}' using RBH only "
                "(RBH produced zero 1:1 orthologs). Pipeline will stop at orthogroup discovery."
            )
        print("[INFO] Keeping RBH result for downstream use (RBH-only mode).")
        return

    print("[INFO] Running OrthoFinder for RBH-vs-OrthoFinder fallback comparison.")

    fallback_dir = outdir / "_orthofinder_fallback"
    try:
        run_orthofinder_fallback(
            query_fasta=query_fasta,
            proteomes_dir=proteomes_dir,
            outdir=fallback_dir,
            coverage=coverage,
            evalue=evalue,
            threads=threads,
            blastp_exe=blastp_exe,
            makeblastdb_exe=makeblastdb_exe,
            orthofinder_exe=orthofinder_exe,
        )
    except Exception as exc:
        raise RuntimeError(
            f"OrthoFinder fallback failed after RBH stage ({rbh_one_to_one}/{total}). "
            f"Pipeline will stop at orthogroup discovery.\nDetails: {exc}"
        ) from exc

    of_rows = read_summary_rows(fallback_dir / "rbh_summary.tsv")
    of_one_to_one = count_one_to_one_orthologs(of_rows)
    of_total = max(1, len(of_rows))
    of_ratio = float(of_one_to_one) / float(of_total)
    print(
        f"[INFO] OrthoFinder 1:1 ortholog recovery for query '{query_id}': "
        f"{of_one_to_one}/{of_total} ({of_ratio:.1%})."
    )

    if max(rbh_one_to_one, of_one_to_one) <= 0:
        raise RuntimeError(
            f"No orthogroup could be defined for query '{query_id}' after RBH and OrthoFinder fallback "
            "(both produced zero 1:1 orthologs). Pipeline will stop at orthogroup discovery."
        )

    if of_one_to_one > rbh_one_to_one:
        for name in ("orthogroup_proteins.fasta", "orthogroup_headers.txt", "rbh_summary.tsv"):
            src = fallback_dir / name
            dst = outdir / name
            if not src.exists():
                raise RuntimeError(f"Missing fallback output file: {src}")
            dst.parent.mkdir(parents=True, exist_ok=True)
            dst.write_bytes(src.read_bytes())
        print(
            f"[INFO] Selected OrthoFinder fallback result for downstream use "
            f"(1:1 orthologs: OrthoFinder={of_one_to_one} > RBH={rbh_one_to_one})."
        )
    elif rbh_one_to_one > of_one_to_one:
        print(
            f"[INFO] Keeping RBH result for downstream use "
            f"(1:1 orthologs: RBH={rbh_one_to_one}, OrthoFinder={of_one_to_one})."
        )
    else:
        print(
            f"[INFO] RBH and OrthoFinder tie on 1:1 ortholog count ({rbh_one_to_one}). "
            "Keeping RBH result for downstream use by deterministic tie-break."
        )


if __name__ == "__main__":
    a = parse_args()
    choose_single_copy_rbh(
        query_fasta=Path(a.query),
        proteomes_dir=Path(a.proteomes),
        outdir=Path(a.outdir),
        coverage=a.coverage,
        evalue=a.evalue,
        threads=a.threads,
        blastp_exe=a.blastp,
        makeblastdb_exe=a.makeblastdb,
        orthofinder_exe=a.orthofinder,
        mode=a.mode,
    )
