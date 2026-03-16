#!/usr/bin/env python3
import argparse
import shutil
import subprocess
from pathlib import Path


MODE_ALIASES = {
    "smartgap": "smart-gap",
    "kpic-smartgap": "kpic-smart-gap",
    "kpic-gappy": "kpic-gappy",
    "kpi-smartgap": "kpi-smart-gap",
    "kpi-gappy": "kpi-gappy",
}


def normalize_mode(mode):
    lowered = mode.strip().lower()
    return MODE_ALIASES.get(lowered, lowered)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", required=True)
    parser.add_argument("--input-fasta", required=True)
    parser.add_argument("--output-fasta", required=True)
    parser.add_argument("--clipkit-log-out", required=True)
    parser.add_argument("--mode", required=True)
    parser.add_argument("--sequence-type", required=True, choices=["aa", "nt"])
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--codon", action="store_true")
    args = parser.parse_args()

    output_fasta = Path(args.output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    clipkit_log_out = Path(args.clipkit_log_out)
    clipkit_log_out.parent.mkdir(parents=True, exist_ok=True)

    staged_input = output_fasta.parent / f"{output_fasta.stem}.clipkit_input{Path(args.input_fasta).suffix or '.fasta'}"
    shutil.copyfile(args.input_fasta, staged_input)
    preexisting_logs = {path.name for path in output_fasta.parent.glob("*.log")}

    command = [
        args.executable,
        staged_input.name,
        "-o",
        output_fasta.name,
        "-m",
        normalize_mode(args.mode),
        "-s",
        args.sequence_type,
        "-l",
        "-t",
        str(args.threads),
        "-q",
    ]
    if args.codon:
        command.append("-co")

    print("[ClipKIT] Running:", " ".join(command))
    result = subprocess.run(
        command,
        cwd=output_fasta.parent,
        capture_output=True,
        text=True,
    )
    if result.stdout.strip():
        print(result.stdout)
    if result.stderr.strip():
        print(result.stderr)
    if result.returncode != 0:
        raise RuntimeError(f"ClipKIT failed with exit code {result.returncode}")
    if not output_fasta.exists():
        raise FileNotFoundError(f"ClipKIT did not create expected output {output_fasta}")

    new_logs = [
        path for path in output_fasta.parent.glob("*.log") if path.name not in preexisting_logs
    ]
    if not new_logs:
        raise FileNotFoundError(
            f"ClipKIT did not create a trimming log alongside {output_fasta}"
        )
    if len(new_logs) > 1:
        new_logs.sort()
    shutil.move(str(new_logs[0]), clipkit_log_out)
    staged_input.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
