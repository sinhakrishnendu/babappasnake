#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path
from typing import Any


def run_command(cmd: list[str]) -> tuple[int, str, str]:
    try:
        res = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return int(res.returncode), res.stdout, res.stderr
    except Exception as exc:  # pragma: no cover - defensive
        return 127, "", f"{type(exc).__name__}: {exc}"


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def write_json(path: Path, payload: dict[str, Any]) -> None:
    write_text(path, json.dumps(payload, indent=2) + "\n")


def _count_breakpoints_from_value(value: Any) -> int | None:
    if value is None:
        return None
    if isinstance(value, int):
        return max(0, int(value))
    if isinstance(value, float):
        if value.is_integer():
            return max(0, int(value))
        return None
    if isinstance(value, (list, tuple)):
        if not value:
            return 0
        numeric = [x for x in value if isinstance(x, (int, float))]
        if len(numeric) == len(value):
            return len(value)
        if all(isinstance(x, dict) for x in value):
            return len(value)
        return None
    if isinstance(value, dict):
        for key in ("count", "n", "n_breakpoints"):
            v = value.get(key)
            if isinstance(v, (int, float)):
                return max(0, int(v))
        if value and all(isinstance(k, (str, int)) for k in value.keys()):
            return len(value)
    return None


def infer_breakpoint_count(payload: Any) -> tuple[int | None, str]:
    candidates: list[int] = []

    def walk(node: Any) -> None:
        if isinstance(node, dict):
            for key, value in node.items():
                key_l = str(key).strip().lower()
                if "breakpoint" in key_l:
                    parsed = _count_breakpoints_from_value(value)
                    if parsed is not None:
                        candidates.append(parsed)
                walk(value)
        elif isinstance(node, list):
            for item in node:
                walk(item)

    walk(payload)
    if not candidates:
        return None, "parse_unavailable"
    return max(candidates), "ok"


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--cds-aln", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--hyphy", default="hyphy")
    p.add_argument("--mode", choices=["none", "gard", "auto"], default="none")
    p.add_argument("--gard-mode", choices=["Normal", "Faster"], default="Faster")
    p.add_argument("--rate-classes", type=int, default=3)
    p.add_argument("--method", default="")
    p.add_argument("--trim-state", default="")
    p.add_argument("--fail-on-error", action="store_true")
    a = p.parse_args()

    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    summary_path = outdir / "gard_summary.json"
    gard_json = outdir / "gard.json"
    stdout_log = outdir / "gard.stdout.txt"
    stderr_log = outdir / "gard.stderr.txt"

    mode_effective = "gard" if a.mode == "auto" else a.mode
    if mode_effective == "none":
        write_json(
            gard_json,
            {
                "status": "not_run",
                "analysis": "gard",
                "note": "Recombination screening mode is 'none'.",
            },
        )
        write_text(stdout_log, "")
        write_text(stderr_log, "")
        write_json(
            summary_path,
            {
                "analysis": "gard",
                "method": a.method,
                "trim_state": a.trim_state,
                "requested_mode": a.mode,
                "effective_mode": mode_effective,
                "status": "not_run",
                "completed": False,
                "breakpoints_detected": False,
                "n_breakpoints": None,
                "parse_status": "not_applicable",
                "primary_json": str(gard_json),
                "stdout_log": str(stdout_log),
                "stderr_log": str(stderr_log),
                "note": "GARD screening not requested.",
            },
        )
        return

    base_cmd = [
        a.hyphy,
        "gard",
        "--alignment",
        a.cds_aln,
        "--output",
        str(gard_json),
    ]
    opt_cmd = base_cmd + [
        "--mode",
        a.gard_mode,
        "--rate-classes",
        str(int(a.rate_classes)),
    ]

    rc, stdout, stderr = run_command(opt_cmd)
    used_fallback = False
    command_used = list(opt_cmd)
    text_for_detect = f"{stdout}\n{stderr}".lower()
    if rc != 0 and ("unknown option" in text_for_detect or "unrecognized option" in text_for_detect):
        fallback_rc, fallback_stdout, fallback_stderr = run_command(base_cmd)
        if fallback_rc == 0:
            rc, stdout, stderr = fallback_rc, fallback_stdout, fallback_stderr
            used_fallback = True
            command_used = list(base_cmd)

    write_text(stdout_log, stdout)
    write_text(stderr_log, stderr)

    if rc != 0:
        write_json(
            gard_json,
            {
                "status": "failed",
                "analysis": "gard",
                "return_code": rc,
                "command": command_used,
                "stdout": stdout,
                "stderr": stderr,
                "fallback_used": used_fallback,
            },
        )
        write_json(
            summary_path,
            {
                "analysis": "gard",
                "method": a.method,
                "trim_state": a.trim_state,
                "requested_mode": a.mode,
                "effective_mode": mode_effective,
                "status": "failed",
                "completed": False,
                "breakpoints_detected": False,
                "n_breakpoints": None,
                "parse_status": "not_applicable",
                "primary_json": str(gard_json),
                "stdout_log": str(stdout_log),
                "stderr_log": str(stderr_log),
                "return_code": rc,
                "command": command_used,
                "fallback_used": used_fallback,
                "note": "HyPhy GARD execution failed.",
            },
        )
        if a.fail_on_error:
            raise RuntimeError(
                f"HyPhy GARD failed with return code {rc}. See {stderr_log} and {gard_json}."
            )
        return

    payload: dict[str, Any] = {}
    parse_status = "parse_unavailable"
    n_breakpoints = None
    try:
        payload = json.loads(gard_json.read_text(encoding="utf-8"))
        n_breakpoints, parse_status = infer_breakpoint_count(payload)
    except Exception:
        parse_status = "parse_unavailable"
        n_breakpoints = None

    breakpoints_detected = bool(n_breakpoints and int(n_breakpoints) > 0)
    note = (
        "GARD completed. Downstream default workflow remains full-length unless fragment-aware routing is explicitly implemented."
    )
    if parse_status != "ok":
        note = (
            "GARD completed, but breakpoint count could not be parsed robustly from JSON schema; recorded as parse_unavailable."
        )

    write_json(
        summary_path,
        {
            "analysis": "gard",
            "method": a.method,
            "trim_state": a.trim_state,
            "requested_mode": a.mode,
            "effective_mode": mode_effective,
            "status": "ok",
            "completed": True,
            "breakpoints_detected": breakpoints_detected,
            "n_breakpoints": n_breakpoints,
            "parse_status": parse_status,
            "primary_json": str(gard_json),
            "stdout_log": str(stdout_log),
            "stderr_log": str(stderr_log),
            "return_code": 0,
            "command": command_used,
            "fallback_used": used_fallback,
            "note": note,
        },
    )


if __name__ == "__main__":
    main()
