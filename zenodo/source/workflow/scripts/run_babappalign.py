#!/usr/bin/env python3
import argparse
import os
import shutil
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", required=True)
    parser.add_argument("--model-path", required=True)
    parser.add_argument("--device", required=True)
    parser.add_argument("--gap-open", required=True, type=float)
    parser.add_argument("--gap-extend", required=True, type=float)
    parser.add_argument("--input-fasta", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--staged-input-out", required=True)
    parser.add_argument("--raw-cds-out", required=True)
    parser.add_argument("--raw-proteins-out", required=True)
    parser.add_argument("--aligned-cds-out", required=True)
    parser.add_argument("--aligned-proteins-out", required=True)
    args = parser.parse_args()

    if not args.model_path:
        raise ValueError(
            "BABAPPAlign requires a scoring model. Set alignment.model_path in config/config.yaml."
        )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    staged_input = Path(args.staged_input_out)
    shutil.copyfile(args.input_fasta, staged_input)
    model_path = Path(args.model_path).expanduser().resolve()
    cache_root = model_path.parents[2] if len(model_path.parents) >= 3 else model_path.parent

    command = [
        args.executable,
        str(staged_input),
        "--model",
        args.model_path,
        "--mode",
        "codon",
        "--device",
        args.device,
        "--gap-open",
        str(args.gap_open),
        "--gap-extend",
        str(args.gap_extend),
    ]
    print("[BABAPPAlign] Running:", " ".join(command))
    env = os.environ.copy()
    env.setdefault("PYTHONNOUSERSITE", "1")
    env.setdefault("HF_HUB_OFFLINE", "1")
    env.setdefault("TRANSFORMERS_OFFLINE", "1")
    env.setdefault("HF_HUB_DISABLE_TELEMETRY", "1")
    env.setdefault("XDG_CACHE_HOME", str(cache_root))
    env.setdefault("HF_HOME", str(cache_root / "huggingface"))
    env.setdefault("TRANSFORMERS_CACHE", str(cache_root / "huggingface" / "hub"))
    result = subprocess.run(command, capture_output=True, text=True, env=env)
    if result.stdout.strip():
        print(result.stdout)
    if result.stderr.strip():
        print(result.stderr)
    if result.returncode != 0:
        raise RuntimeError(f"BABAPPAlign failed with exit code {result.returncode}")

    protein_candidates = [
        staged_input.with_suffix(staged_input.suffix + ".protein.aln.fasta"),
        staged_input.with_suffix(".protein.aln.fasta"),
    ]
    codon_candidates = [
        staged_input.with_suffix(staged_input.suffix + ".codon.aln.fasta"),
        staged_input.with_suffix(".codon.aln.fasta"),
    ]
    raw_protein = next((path for path in protein_candidates if path.exists()), None)
    raw_codon = next((path for path in codon_candidates if path.exists()), None)
    if raw_protein is None or raw_codon is None:
        raise FileNotFoundError(
            "BABAPPAlign finished without producing the expected codon and protein alignments."
        )

    for source, destination in (
        (raw_protein, Path(args.raw_proteins_out)),
        (raw_codon, Path(args.raw_cds_out)),
        (raw_protein, Path(args.aligned_proteins_out)),
        (raw_codon, Path(args.aligned_cds_out)),
    ):
        if Path(source).resolve() != destination.resolve():
            shutil.copyfile(source, destination)


if __name__ == "__main__":
    main()
