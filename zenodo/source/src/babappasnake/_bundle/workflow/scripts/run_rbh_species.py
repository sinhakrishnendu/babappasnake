#!/usr/bin/env python3
import argparse
import csv
import subprocess
from pathlib import Path

from common import guess_taxon_from_identifier, read_fasta_records


SEARCH_COLUMNS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore"]


def run_command(command, log_prefix):
    result = subprocess.run(command, capture_output=True, text=True)
    if result.stdout.strip():
        print(f"[{log_prefix} stdout]\n{result.stdout}")
    if result.stderr.strip():
        print(f"[{log_prefix} stderr]\n{result.stderr}")
    if result.returncode != 0:
        raise RuntimeError(f"{log_prefix} failed with exit code {result.returncode}")


def ensure_blast_db(makeblastdb, fasta, prefix):
    if Path(prefix + ".pin").exists():
        return
    run_command(
        [makeblastdb, "-in", str(fasta), "-dbtype", "prot", "-out", prefix],
        "makeblastdb",
    )


def ensure_diamond_db(diamond, fasta, prefix):
    if Path(prefix + ".dmnd").exists():
        return
    run_command([diamond, "makedb", "--in", str(fasta), "--db", prefix], "diamond-makedb")


def run_search(args, query_fasta, db_prefix, output_path):
    if args.search_tool == "blast":
        command = [
            args.blastp,
            "-query",
            str(query_fasta),
            "-db",
            db_prefix,
            "-out",
            str(output_path),
            "-evalue",
            str(args.evalue),
            "-outfmt",
            "6 qseqid sseqid pident length qlen slen evalue bitscore",
            "-max_target_seqs",
            str(args.max_target_seqs),
            "-num_threads",
            str(args.threads),
            "-seg",
            "yes",
            "-soft_masking",
            "true",
            "-use_sw_tback",
        ]
    else:
        command = [
            args.diamond,
            "blastp",
            "--query",
            str(query_fasta),
            "--db",
            db_prefix,
            "--out",
            str(output_path),
            "--evalue",
            str(args.evalue),
            "--max-target-seqs",
            str(args.max_target_seqs),
            "--threads",
            str(args.threads),
            "--outfmt",
            "6",
            *SEARCH_COLUMNS,
        ]
        if args.diamond_mode:
            command.append(f"--{args.diamond_mode}")
    run_command(command, f"{args.search_tool}-search")


def read_search_table(path):
    hits = []
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t", fieldnames=SEARCH_COLUMNS)
        for row in reader:
            row["pident"] = float(row["pident"])
            row["length"] = int(float(row["length"]))
            row["qlen"] = int(float(row["qlen"]))
            row["slen"] = int(float(row["slen"]))
            row["evalue"] = float(row["evalue"])
            row["bitscore"] = float(row["bitscore"])
            hits.append(row)
    return hits


def best_hits_by_query(hits):
    grouped = {}
    for row in hits:
        grouped.setdefault(row["qseqid"], []).append(row)
    for value in grouped.values():
        value.sort(key=lambda row: (-row["bitscore"], row["evalue"], row["sseqid"]))
    return grouped


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--query-fasta", required=True)
    parser.add_argument("--proteome-fasta", required=True)
    parser.add_argument("--species-id", required=True)
    parser.add_argument("--taxon-name", required=True)
    parser.add_argument("--search-tool", choices=["blast", "diamond"], required=True)
    parser.add_argument("--blastp", default="blastp")
    parser.add_argument("--makeblastdb", default="makeblastdb")
    parser.add_argument("--diamond", default="diamond")
    parser.add_argument("--diamond-mode", default="")
    parser.add_argument("--evalue", required=True, type=float)
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--max-target-seqs", required=True, type=int)
    parser.add_argument("--query-db-prefix", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--forward-out", required=True)
    parser.add_argument("--reverse-out", required=True)
    parser.add_argument("--candidate-out", required=True)
    args = parser.parse_args()

    query_records = read_fasta_records(args.query_fasta)
    if len(query_records) != 1:
        raise ValueError("run_rbh_species.py expects a single-sequence query FASTA")
    query_id = query_records[0]["id"]

    workdir = Path(args.workdir)
    db_dir = workdir / "db"
    db_dir.mkdir(parents=True, exist_ok=True)

    if args.search_tool == "blast":
        ensure_blast_db(args.makeblastdb, args.query_fasta, args.query_db_prefix)
        proteome_db_prefix = str(db_dir / "proteome")
        ensure_blast_db(args.makeblastdb, args.proteome_fasta, proteome_db_prefix)
    else:
        ensure_diamond_db(args.diamond, args.query_fasta, args.query_db_prefix)
        proteome_db_prefix = str(db_dir / "proteome")
        ensure_diamond_db(args.diamond, args.proteome_fasta, proteome_db_prefix)

    forward_out = Path(args.forward_out)
    reverse_out = Path(args.reverse_out)
    forward_out.parent.mkdir(parents=True, exist_ok=True)

    run_search(args, args.query_fasta, proteome_db_prefix, forward_out)
    run_search(args, args.proteome_fasta, args.query_db_prefix, reverse_out)

    forward_hits = read_search_table(forward_out)
    reverse_hits = read_search_table(reverse_out)
    reverse_best = {
        query: rows[0]
        for query, rows in best_hits_by_query(reverse_hits).items()
    }

    query_forward_hits = best_hits_by_query(forward_hits).get(query_id, [])
    candidate_rows = []
    for rank, row in enumerate(query_forward_hits, start=1):
        reverse_best_hit = reverse_best.get(row["sseqid"])
        if not reverse_best_hit or reverse_best_hit["sseqid"] != query_id:
            continue
        candidate_taxon = guess_taxon_from_identifier(row["sseqid"]) or args.taxon_name
        candidate_rows.append(
            {
                "species_id": args.species_id,
                "taxon_name": candidate_taxon,
                "query_id": row["qseqid"],
                "subject_id": row["sseqid"],
                "rank_by_bitscore": rank,
                "bitscore": row["bitscore"],
                "evalue": row["evalue"],
                "pident": row["pident"],
                "alignment_length": row["length"],
                "query_length": row["qlen"],
                "subject_length": row["slen"],
                "query_coverage": row["length"] / row["qlen"],
                "subject_coverage": row["length"] / row["slen"],
                "reverse_bitscore": reverse_best_hit["bitscore"],
                "proteome_fasta": str(Path(args.proteome_fasta).resolve()),
            }
        )

    candidate_path = Path(args.candidate_out)
    candidate_path.parent.mkdir(parents=True, exist_ok=True)
    with open(candidate_path, "w", newline="") as handle:
        fieldnames = [
            "species_id",
            "taxon_name",
            "query_id",
            "subject_id",
            "rank_by_bitscore",
            "bitscore",
            "evalue",
            "pident",
            "alignment_length",
            "query_length",
            "subject_length",
            "query_coverage",
            "subject_coverage",
            "reverse_bitscore",
            "proteome_fasta",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(candidate_rows)

    print(
        f"[RBH] {args.species_id}: retained {len(candidate_rows)} reciprocal candidate(s) "
        f"from {len(query_forward_hits)} forward hit(s)"
    )


if __name__ == "__main__":
    main()
