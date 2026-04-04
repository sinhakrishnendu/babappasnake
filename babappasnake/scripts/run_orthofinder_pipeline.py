#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from babappasnake.scripts.run_rbh_pipeline import discover_proteome_fastas


QUERY_SPECIES_LABEL = "babappasnake_query"
BLAST_SUBJECT_SEP = "||"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--query", required=True)
    p.add_argument("--proteomes", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--coverage", type=float, default=0.70)
    p.add_argument("--evalue", type=float, default=1e-5)
    p.add_argument("--blastp", default="blastp")
    p.add_argument("--makeblastdb", default="makeblastdb")
    p.add_argument("--orthofinder", default="orthofinder")
    return p.parse_args()


def run_cmd(cmd: list[str]) -> None:
    result = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed:\n{' '.join(cmd)}\n\nSTDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}"
        )


def parse_members(cell: str) -> list[str]:
    value = str(cell or "").strip()
    if not value:
        return []
    return [token.strip() for token in value.split(",") if token.strip()]


def find_latest_orthofinder_results(input_dir: Path) -> Path:
    root = input_dir / "OrthoFinder"
    if not root.exists():
        raise RuntimeError(f"OrthoFinder output directory not found: {root}")
    candidates = [p for p in root.glob("Results_*") if p.is_dir()]
    if not candidates:
        raise RuntimeError(f"No OrthoFinder Results_* directory found under: {root}")
    return max(candidates, key=lambda p: p.stat().st_mtime)


def load_orthogroups_from_tsv(
    orthogroups_tsv: Path,
    query_species: str,
) -> tuple[list[str], Dict[str, Dict[str, list[str]]]]:
    if not orthogroups_tsv.exists():
        raise RuntimeError(f"Orthogroups table not found: {orthogroups_tsv}")

    with open(orthogroups_tsv, encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        if "Orthogroup" not in fieldnames:
            raise RuntimeError(f"Orthogroups table missing required 'Orthogroup' column: {orthogroups_tsv}")
        if query_species not in fieldnames:
            raise RuntimeError(
                f"Orthogroups table missing query species column '{query_species}'. "
                f"Columns: {', '.join(fieldnames)}"
            )

        species_columns = [name for name in fieldnames if name not in {"Orthogroup", query_species}]
        groups: Dict[str, Dict[str, list[str]]] = {}
        for row in reader:
            orthogroup_id = str(row.get("Orthogroup", "")).strip()
            if not orthogroup_id:
                continue
            members_by_species: Dict[str, list[str]] = {}
            for species in species_columns:
                members_by_species[species] = parse_members(row.get(species, ""))
            groups[orthogroup_id] = members_by_species

    if not groups:
        raise RuntimeError(f"No orthogroups were found in OrthoFinder table: {orthogroups_tsv}")
    return species_columns, groups


def build_sequence_index(proteome_fastas: list[Path]) -> tuple[Dict[str, SeqRecord], set[str]]:
    index: Dict[str, SeqRecord] = {}
    duplicates: set[str] = set()
    for fasta in proteome_fastas:
        for rec in SeqIO.parse(str(fasta), "fasta"):
            if rec.id in index:
                duplicates.add(rec.id)
                continue
            index[rec.id] = rec
    return index, duplicates


def build_orthogroup_subject_fasta(
    outdir: Path,
    species_columns: list[str],
    orthogroups: Dict[str, Dict[str, list[str]]],
    seq_index: Dict[str, SeqRecord],
    duplicates: set[str],
) -> Path:
    search_dir = outdir / "orthofinder_query_search"
    search_dir.mkdir(parents=True, exist_ok=True)
    subject_fasta = search_dir / "orthogroups_subjects.fasta"

    subject_records: list[SeqRecord] = []
    seen_subject_ids: set[str] = set()
    for orthogroup_id in sorted(orthogroups):
        members_by_species = orthogroups[orthogroup_id]
        for species in species_columns:
            for member in members_by_species.get(species, []):
                if member in duplicates:
                    raise RuntimeError(
                        f"Orthogroup {orthogroup_id}: ambiguous duplicated sequence ID across proteomes: {member}"
                    )
                rec = seq_index.get(member)
                if rec is None:
                    raise RuntimeError(
                        f"Orthogroup {orthogroup_id}: could not resolve sequence ID in proteomes: {member}"
                    )
                subject_id = f"{orthogroup_id}{BLAST_SUBJECT_SEP}{species}{BLAST_SUBJECT_SEP}{member}"
                if subject_id in seen_subject_ids:
                    continue
                seen_subject_ids.add(subject_id)
                subject_records.append(SeqRecord(rec.seq, id=subject_id, description=""))

    if not subject_records:
        raise RuntimeError(
            "OrthoFinder produced no orthogroup partner sequences for BLAST-based query mapping."
        )

    SeqIO.write(subject_records, subject_fasta, "fasta")
    return subject_fasta


def make_blast_db(fasta: Path, makeblastdb_exe: str) -> str:
    prefix = str(fasta.with_suffix(""))
    if not Path(prefix + ".pin").exists():
        run_cmd([makeblastdb_exe, "-in", str(fasta), "-dbtype", "prot", "-out", prefix])
    return prefix


def run_query_blast(
    query_fasta: Path,
    db_prefix: str,
    out_tsv: Path,
    threads: int,
    evalue: float,
    blastp_exe: str,
) -> None:
    run_cmd(
        [
            blastp_exe,
            "-query",
            str(query_fasta),
            "-db",
            db_prefix,
            "-out",
            str(out_tsv),
            "-evalue",
            str(evalue),
            "-outfmt",
            "6 qseqid sseqid pident length qlen slen evalue bitscore",
            "-num_threads",
            str(max(1, int(threads))),
            "-seg",
            "yes",
            "-soft_masking",
            "true",
            "-use_sw_tback",
        ]
    )


def select_best_orthogroup_from_blast_tsv(blast_tsv: Path, min_coverage: float) -> str:
    if not blast_tsv.exists():
        raise RuntimeError(f"BLAST result file not found: {blast_tsv}")

    grouped_scores: Dict[str, Dict[str, float]] = {}
    grouped_best: Dict[str, float] = {}
    with open(blast_tsv, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            sseqid = parts[1]
            length = int(parts[3])
            qlen = int(parts[4])
            bitscore = float(parts[7])
            if qlen <= 0:
                continue
            qcov = length / qlen
            if qcov < min_coverage:
                continue
            tokens = sseqid.split(BLAST_SUBJECT_SEP, 2)
            if len(tokens) != 3:
                continue
            orthogroup_id, species, _member = tokens
            species_scores = grouped_scores.setdefault(orthogroup_id, {})
            prev = species_scores.get(species)
            if prev is None or bitscore > prev:
                species_scores[species] = bitscore
            grouped_best[orthogroup_id] = max(grouped_best.get(orthogroup_id, 0.0), bitscore)

    if not grouped_scores:
        raise RuntimeError(
            f"No BLAST hits passed coverage >= {min_coverage:.2f} when mapping query to OrthoFinder orthogroups."
        )

    ranked: list[tuple[int, float, float, str]] = []
    for orthogroup_id, species_scores in grouped_scores.items():
        species_hit_count = len(species_scores)
        total_bitscore = sum(species_scores.values())
        top_bitscore = grouped_best.get(orthogroup_id, 0.0)
        ranked.append((species_hit_count, total_bitscore, top_bitscore, orthogroup_id))

    ranked.sort(key=lambda item: (-item[0], -item[1], -item[2], item[3]))
    return ranked[0][3]


def prepare_orthofinder_input(
    outdir: Path,
    proteome_fastas: list[Path],
    query_record: SeqRecord,
) -> Path:
    input_dir = outdir / "orthofinder_input"
    if input_dir.exists():
        shutil.rmtree(input_dir)
    input_dir.mkdir(parents=True, exist_ok=True)

    for fasta in proteome_fastas:
        dst = input_dir / fasta.name
        try:
            dst.symlink_to(fasta.resolve())
        except OSError:
            shutil.copy2(fasta, dst)

    query_input = input_dir / f"{QUERY_SPECIES_LABEL}.fasta"
    SeqIO.write([query_record], query_input, "fasta")
    return input_dir


def choose_orthogroup_with_orthofinder(
    query_fasta: Path,
    proteomes_dir: Path,
    outdir: Path,
    threads: int,
    coverage: float,
    evalue: float,
    blastp_exe: str,
    makeblastdb_exe: str,
    orthofinder_exe: str,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    query_records = list(SeqIO.parse(str(query_fasta), "fasta"))
    if len(query_records) != 1:
        raise ValueError("OrthoFinder mode expects exactly one query protein in query FASTA.")
    query_record = query_records[0]
    query_id = query_record.id

    proteome_fastas = discover_proteome_fastas(proteomes_dir)
    seq_index, duplicates = build_sequence_index(proteome_fastas)

    input_dir = prepare_orthofinder_input(outdir, proteome_fastas, query_record)
    run_cmd(
        [
            orthofinder_exe,
            "-f",
            str(input_dir),
            "-t",
            str(max(1, int(threads))),
        ]
    )
    results_dir = find_latest_orthofinder_results(input_dir)
    orthogroups_tsv = results_dir / "Orthogroups" / "Orthogroups.tsv"
    species_columns, orthogroups = load_orthogroups_from_tsv(
        orthogroups_tsv=orthogroups_tsv,
        query_species=QUERY_SPECIES_LABEL,
    )

    subject_fasta = build_orthogroup_subject_fasta(
        outdir=outdir,
        species_columns=species_columns,
        orthogroups=orthogroups,
        seq_index=seq_index,
        duplicates=duplicates,
    )
    db_prefix = make_blast_db(subject_fasta, makeblastdb_exe)
    blast_tsv = subject_fasta.parent / "query_vs_orthogroups.tsv"
    run_query_blast(
        query_fasta=query_fasta,
        db_prefix=db_prefix,
        out_tsv=blast_tsv,
        threads=threads,
        evalue=evalue,
        blastp_exe=blastp_exe,
    )
    orthogroup_id = select_best_orthogroup_from_blast_tsv(blast_tsv, min_coverage=coverage)
    members_by_species = orthogroups.get(orthogroup_id, {})
    if not members_by_species:
        raise RuntimeError(
            f"Selected orthogroup '{orthogroup_id}' from BLAST mapping is missing from {orthogroups_tsv}."
        )
    print(
        f"[INFO] OrthoFinder target orthogroup selected by BLAST for query '{query_id}': {orthogroup_id}",
        file=sys.stderr,
    )

    summary_rows: list[tuple[str, str, str]] = []
    ortholog_records: list[SeqRecord] = [query_record]
    headers: list[str] = [query_id]
    seen_headers = {query_id}

    for species in species_columns:
        members = members_by_species.get(species, [])
        if len(members) != 1:
            if len(members) > 1:
                print(
                    f"[INFO] OrthoFinder species '{species}' has {len(members)} ortholog candidates in "
                    f"orthogroup {orthogroup_id}; excluding from 1:1 set.",
                    file=sys.stderr,
                )
            summary_rows.append((species, query_id, "NA"))
            continue

        member = members[0]
        if member not in seq_index:
            raise RuntimeError(
                f"Orthogroup {orthogroup_id}: could not resolve sequence ID in proteomes: {member}"
            )

        if member in duplicates:
            raise RuntimeError(
                f"Orthogroup {orthogroup_id}: ambiguous duplicated sequence ID across proteomes: {member}"
            )

        if member in seen_headers:
            summary_rows.append((species, query_id, "NA"))
            continue
        ortholog_records.append(seq_index[member])
        headers.append(member)
        seen_headers.add(member)
        summary_rows.append((species, query_id, member))

    partner_count = len(ortholog_records) - 1
    if partner_count <= 0:
        raise RuntimeError(
            f"No orthogroup could be defined for query '{query_id}' using OrthoFinder: "
            f"BLAST-selected orthogroup '{orthogroup_id}' has no partner sequences from the provided proteomes. "
            "Pipeline will stop at orthogroup discovery."
        )

    SeqIO.write(ortholog_records, outdir / "orthogroup_proteins.fasta", "fasta")
    with open(outdir / "orthogroup_headers.txt", "w", encoding="utf-8") as fh:
        for header in headers:
            fh.write(header + "\n")
    with open(outdir / "rbh_summary.tsv", "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["species", "query", "ortholog"])
        writer.writerows(summary_rows)


if __name__ == "__main__":
    a = parse_args()
    choose_orthogroup_with_orthofinder(
        query_fasta=Path(a.query),
        proteomes_dir=Path(a.proteomes),
        outdir=Path(a.outdir),
        threads=int(a.threads),
        coverage=float(a.coverage),
        evalue=float(a.evalue),
        blastp_exe=str(a.blastp),
        makeblastdb_exe=str(a.makeblastdb),
        orthofinder_exe=str(a.orthofinder),
    )
