#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from itertools import combinations
from pathlib import Path

from Bio import SeqIO


def parse_pathways(raw: str) -> list[tuple[str, str]]:
    pathways: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for token in [x.strip() for x in str(raw).split(",") if x.strip()]:
        method = ""
        trim_state = ""
        if ":" in token:
            method, trim_state = token.split(":", 1)
        elif token.endswith("_raw"):
            method, trim_state = token[:-4], "raw"
        elif token.endswith("_clipkit"):
            method, trim_state = token[:-8], "clipkit"
        else:
            method, trim_state = token, "clipkit"
        key = (method.strip(), trim_state.strip())
        if not all(key):
            continue
        if key not in seen:
            seen.add(key)
            pathways.append(key)
    if not pathways:
        raise RuntimeError("No pathways were provided.")
    return pathways


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def parse_bool(value: str) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def detect_significance_key(rows: list[dict[str, str]]) -> str | None:
    preferred = ("significant_BH_0.05", "significant_BH_0.1", "significant_BH")
    available = {k for row in rows for k in row.keys()}
    for key in preferred:
        if key in available:
            return key
    for key in sorted(available):
        if key.startswith("significant_BH_"):
            return key
    return None


def read_fasta_shape(path: Path) -> tuple[int, int]:
    if not path.exists() or path.stat().st_size == 0:
        return 0, 0
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        return 0, 0
    length = len(str(records[0].seq))
    return len(records), int(length)


def read_absrel_tested_branches(path: Path) -> int:
    if not path.exists() or path.stat().st_size == 0:
        return 0
    data = json.loads(path.read_text(encoding="utf-8"))
    fits = data.get("branch attributes", {}).get("0", {})
    if not isinstance(fits, dict):
        return 0
    count = 0
    for _branch, stats in fits.items():
        if not isinstance(stats, dict):
            continue
        is_leaf = stats.get("is leaf", None)
        if is_leaf is True or is_leaf is None:
            count += 1
    return count


def extract_meme_hits(payload: object, threshold: float) -> tuple[int, set[str]]:
    count = 0
    sites: set[str] = set()

    def walk(node: object, inherited_site: str | None = None) -> None:
        nonlocal count
        if isinstance(node, dict):
            site = inherited_site
            for key in (
                "site",
                "Site",
                "codon",
                "Codon",
                "position",
                "Position",
                "site_index",
                "Site Index",
            ):
                if key in node:
                    site = str(node.get(key)).strip()
                    break

            pval = None
            for key in ("p-value", "P-value", "p", "pvalue"):
                if key in node:
                    try:
                        pval = float(node.get(key))
                    except Exception:
                        pval = None
                    break

            if pval is not None and pval <= threshold:
                count += 1
                if site:
                    sites.add(site)

            for value in node.values():
                walk(value, site)

        elif isinstance(node, list):
            for idx, value in enumerate(node, start=1):
                next_site = inherited_site if inherited_site is not None else str(idx)
                walk(value, next_site)

    walk(payload)
    return count, sites


def jaccard(a: set[str], b: set[str]) -> float:
    union = a | b
    if not union:
        return 1.0
    return len(a & b) / len(union)


def classify_signal(
    replication_count: int,
    total_pathways: int,
    method_set: set[str],
    trim_set: set[str],
    total_methods: int,
    total_trim_states: int,
) -> str:
    if replication_count == total_pathways:
        return "highly_robust"

    method_specific = total_methods > 1 and len(method_set) == 1
    trim_specific = total_trim_states > 1 and len(trim_set) == 1
    if method_specific and not trim_specific:
        return "method_sensitive"
    if trim_specific and not method_specific:
        return "trim_sensitive"
    if total_pathways > 0 and (replication_count / total_pathways) >= 0.5:
        return "moderately_robust"
    return "not_reproducible"


def escape_latex(text: str) -> str:
    mapping = {
        "\\": r"\textbackslash{}",
        "_": r"\_",
        "%": r"\%",
        "&": r"\&",
        "#": r"\#",
        "{": r"\{",
        "}": r"\}",
        "$": r"\$",
    }
    return "".join(mapping.get(ch, ch) for ch in str(text))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", required=True)
    p.add_argument("--pathways", required=True)
    p.add_argument("--meme-p", type=float, default=0.05)
    p.add_argument("--matrix-out", required=True)
    p.add_argument("--consensus-out", required=True)
    p.add_argument("--narrative-out", required=True)
    p.add_argument("--comparative-out", required=True)
    p.add_argument("--latex-out", required=True)
    a = p.parse_args()

    outdir = Path(a.outdir)
    pathways = parse_pathways(a.pathways)
    pathway_names = [f"{m}_{t}" for m, t in pathways]
    pathway_to_parts = {f"{m}_{t}": (m, t) for m, t in pathways}
    method_set = {m for m, _ in pathways}
    trim_set = {t for _, t in pathways}

    matrix_rows: list[dict[str, str]] = []
    absrel_hits_by_pathway: dict[str, set[str]] = {}
    branchsite_hits_by_pathway: dict[str, set[str]] = {}
    meme_hits_by_pathway: dict[str, set[str]] = {}

    for method, trim_state in pathways:
        pathway = f"{method}_{trim_state}"
        protein_path = outdir / "alignments" / method / trim_state / "orthogroup_proteins.analysis.fasta"
        cds_path = outdir / "alignments" / method / trim_state / "mapped_orthogroup_cds.analysis.fasta"
        rooted_tree = outdir / "tree" / method / trim_state / "orthogroup.rooted.treefile"
        absrel_tsv = outdir / "hyphy" / method / trim_state / "significant_foregrounds.tsv"
        absrel_meta = outdir / "hyphy" / method / trim_state / "foreground_threshold.json"
        absrel_json = outdir / "hyphy" / method / trim_state / "absrel.json"
        meme_json = outdir / "hyphy" / method / trim_state / "meme.json"
        branchsite_tsv = outdir / "branchsite" / method / trim_state / "branchsite_results.tsv"
        asr_done = outdir / "asr" / method / trim_state / "asr_done.json"
        hyphy_done = outdir / "hyphy" / method / trim_state / "hyphy_done.json"

        n_seq, protein_len = read_fasta_shape(protein_path)
        _, cds_len = read_fasta_shape(cds_path)

        absrel_rows = read_tsv(absrel_tsv)
        absrel_sig = sorted({str(row.get("foreground_branch", "")).strip() for row in absrel_rows if str(row.get("foreground_branch", "")).strip()})
        absrel_hits_by_pathway[pathway] = set(absrel_sig)

        branchsite_rows = read_tsv(branchsite_tsv)
        signif_key = detect_significance_key(branchsite_rows)
        branchsite_sig = sorted(
            {
                str(row.get("foreground_branch", "")).strip()
                for row in branchsite_rows
                if str(row.get("foreground_branch", "")).strip()
                and (parse_bool(str(row.get(signif_key, ""))) if signif_key else False)
            }
        )
        branchsite_hits_by_pathway[pathway] = set(branchsite_sig)

        meme_hit_count = 0
        meme_sites: set[str] = set()
        if meme_json.exists() and meme_json.stat().st_size > 0:
            try:
                meme_payload = json.loads(meme_json.read_text(encoding="utf-8"))
                meme_hit_count, meme_sites = extract_meme_hits(meme_payload, float(a.meme_p))
            except Exception:
                meme_hit_count, meme_sites = 0, set()
        meme_hits_by_pathway[pathway] = set(meme_sites)

        absrel_cutoff = "NA"
        if absrel_meta.exists() and absrel_meta.stat().st_size > 0:
            try:
                absrel_cutoff = str(json.loads(absrel_meta.read_text(encoding="utf-8")).get("selected_pcut", "NA"))
            except Exception:
                absrel_cutoff = "NA"

        missing = []
        required_map = {
            "protein_alignment": protein_path,
            "codon_alignment": cds_path,
            "rooted_tree": rooted_tree,
            "absrel_tsv": absrel_tsv,
            "branchsite_tsv": branchsite_tsv,
            "asr_done": asr_done,
            "hyphy_done": hyphy_done,
        }
        for key, path in required_map.items():
            if not path.exists() or path.stat().st_size == 0:
                missing.append(key)

        hyphy_failed = False
        if hyphy_done.exists() and hyphy_done.stat().st_size > 0:
            try:
                hd = json.loads(hyphy_done.read_text(encoding="utf-8"))
                hyphy_failed = (
                    str(hd.get("absrel_status", "ok")).lower() != "ok"
                    or str(hd.get("meme_status", "ok")).lower() != "ok"
                )
            except Exception:
                hyphy_failed = False

        if missing:
            status = f"incomplete_missing:{','.join(missing)}"
        elif hyphy_failed:
            status = "hyphy_failed"
        else:
            status = "ok"

        matrix_rows.append(
            {
                "method": method,
                "trim_state": trim_state,
                "pathway": pathway,
                "n_sequences": str(n_seq),
                "protein_alignment_length": str(protein_len),
                "codon_alignment_length": str(cds_len),
                "rooted_tree_generated": "yes" if rooted_tree.exists() and rooted_tree.stat().st_size > 0 else "no",
                "n_absrel_tested_branches": str(read_absrel_tested_branches(absrel_json)),
                "n_absrel_significant": str(len(absrel_sig)),
                "absrel_significant_branches": ";".join(absrel_sig),
                "n_meme_sites_p_lt_threshold": str(meme_hit_count),
                "n_branchsite_foregrounds_tested": str(len(branchsite_rows)),
                "n_branchsite_significant_after_bh": str(len(branchsite_sig)),
                "branchsite_significant_foregrounds": ";".join(branchsite_sig),
                "asr_completed": "yes" if asr_done.exists() and asr_done.stat().st_size > 0 else "no",
                "status": status,
                "absrel_threshold": absrel_cutoff,
                "meme_threshold": str(a.meme_p),
            }
        )

    matrix_fields = [
        "method",
        "trim_state",
        "pathway",
        "n_sequences",
        "protein_alignment_length",
        "codon_alignment_length",
        "rooted_tree_generated",
        "n_absrel_tested_branches",
        "n_absrel_significant",
        "absrel_significant_branches",
        "n_meme_sites_p_lt_threshold",
        "n_branchsite_foregrounds_tested",
        "n_branchsite_significant_after_bh",
        "branchsite_significant_foregrounds",
        "asr_completed",
        "status",
        "absrel_threshold",
        "meme_threshold",
    ]
    write_tsv(Path(a.matrix_out), matrix_fields, matrix_rows)

    signal_support: dict[tuple[str, str], set[str]] = defaultdict(set)
    for pathway, hits in absrel_hits_by_pathway.items():
        for branch in hits:
            signal_support[("aBSREL_branch", branch)].add(pathway)
    for pathway, hits in branchsite_hits_by_pathway.items():
        for branch in hits:
            signal_support[("branchsite_foreground", branch)].add(pathway)
    for pathway, hits in meme_hits_by_pathway.items():
        for site in hits:
            signal_support[("MEME_site", site)].add(pathway)

    total_pathways = len(pathway_names)
    consensus_rows: list[dict[str, str]] = []
    for (signal_source, signal_id), supporting in sorted(signal_support.items(), key=lambda x: (x[0][0], x[0][1])):
        present_methods = {pathway_to_parts[p][0] for p in supporting}
        present_trim_states = {pathway_to_parts[p][1] for p in supporting}
        label = classify_signal(
            replication_count=len(supporting),
            total_pathways=total_pathways,
            method_set=present_methods,
            trim_set=present_trim_states,
            total_methods=len(method_set),
            total_trim_states=len(trim_set),
        )
        consensus_rows.append(
            {
                "signal_source": signal_source,
                "signal_id": signal_id,
                "replication_count": str(len(supporting)),
                "replication_fraction": f"{(len(supporting) / total_pathways):.3f}" if total_pathways else "0.000",
                "methods_present": ",".join(sorted(present_methods)),
                "trim_states_present": ",".join(sorted(present_trim_states)),
                "pathways_present": ",".join(sorted(supporting)),
                "interpretation": label,
            }
        )

    consensus_fields = [
        "signal_source",
        "signal_id",
        "replication_count",
        "replication_fraction",
        "methods_present",
        "trim_states_present",
        "pathways_present",
        "interpretation",
    ]
    write_tsv(Path(a.consensus_out), consensus_fields, consensus_rows)

    # Sensitivity metrics by method and trimming dimensions.
    lines: list[str] = []
    lines.append("BABAPPASNAKE ROBUSTNESS COMPARATIVE SUMMARY")
    lines.append("=" * 46)
    lines.append("")
    lines.append(f"Pathways evaluated ({total_pathways}): {', '.join(pathway_names)}")
    lines.append("")

    def _avg_jaccard(grouped: list[set[str]]) -> float:
        if len(grouped) < 2:
            return 1.0
        vals = [jaccard(a_set, b_set) for a_set, b_set in combinations(grouped, 2)]
        return sum(vals) / len(vals) if vals else 1.0

    lines.append("Alignment-method sensitivity")
    lines.append("-" * 28)
    for trim_state in sorted(trim_set):
        absrel_sets = [absrel_hits_by_pathway.get(f"{method}_{trim_state}", set()) for method in sorted(method_set) if f"{method}_{trim_state}" in absrel_hits_by_pathway]
        branch_sets = [branchsite_hits_by_pathway.get(f"{method}_{trim_state}", set()) for method in sorted(method_set) if f"{method}_{trim_state}" in branchsite_hits_by_pathway]
        lines.append(
            f"- trim_state={trim_state}: avg Jaccard aBSREL={_avg_jaccard(absrel_sets):.3f}, "
            f"branch-site={_avg_jaccard(branch_sets):.3f}"
        )
    lines.append("")

    lines.append("Trim sensitivity")
    lines.append("-" * 16)
    for method in sorted(method_set):
        absrel_sets = [absrel_hits_by_pathway.get(f"{method}_{trim_state}", set()) for trim_state in sorted(trim_set) if f"{method}_{trim_state}" in absrel_hits_by_pathway]
        branch_sets = [branchsite_hits_by_pathway.get(f"{method}_{trim_state}", set()) for trim_state in sorted(trim_set) if f"{method}_{trim_state}" in branchsite_hits_by_pathway]
        lines.append(
            f"- method={method}: avg Jaccard aBSREL={_avg_jaccard(absrel_sets):.3f}, "
            f"branch-site={_avg_jaccard(branch_sets):.3f}"
        )

    label_counts: dict[str, int] = defaultdict(int)
    for row in consensus_rows:
        label_counts[row["interpretation"]] += 1

    lines.append("")
    lines.append("Signal reproducibility labels")
    lines.append("-" * 29)
    for label in ("highly_robust", "moderately_robust", "method_sensitive", "trim_sensitive", "not_reproducible"):
        lines.append(f"- {label}: {label_counts.get(label, 0)}")

    Path(a.comparative_out).write_text("\n".join(lines) + "\n", encoding="utf-8")

    # Narrative synthesis.
    narrative: list[str] = []
    narrative.append("BABAPPASNAKE ROBUSTNESS NARRATIVE")
    narrative.append("=" * 33)
    narrative.append("")

    for source_name, source_map in (
        ("aBSREL branches", absrel_hits_by_pathway),
        ("branch-site foregrounds", branchsite_hits_by_pathway),
        ("MEME sites", meme_hits_by_pathway),
    ):
        raw_union = set()
        clip_union = set()
        for method in method_set:
            raw_union |= source_map.get(f"{method}_raw", set())
            clip_union |= source_map.get(f"{method}_clipkit", set())

        stable = raw_union & clip_union
        trim_only = clip_union - raw_union
        raw_only = raw_union - clip_union

        narrative.append(f"{source_name}:")
        narrative.append(f"- recovered regardless of trimming: {len(stable)}")
        narrative.append(f"- recovered only after trimming: {len(trim_only)}")
        narrative.append(f"- recovered only without trimming: {len(raw_only)}")

    highly = [row for row in consensus_rows if row["interpretation"] == "highly_robust"]
    method_sensitive = [row for row in consensus_rows if row["interpretation"] == "method_sensitive"]
    trim_sensitive = [row for row in consensus_rows if row["interpretation"] == "trim_sensitive"]

    narrative.append("")
    narrative.append("Overall interpretation")
    narrative.append(f"- highly robust signals across alignment + trimming: {len(highly)}")
    narrative.append(f"- method-sensitive signals: {len(method_sensitive)}")
    narrative.append(f"- trim-sensitive signals: {len(trim_sensitive)}")
    narrative.append(
        "- conclusions are most stable when a signal is repeatedly recovered across both trim states and multiple alignment engines."
    )

    Path(a.narrative_out).write_text("\n".join(narrative) + "\n", encoding="utf-8")

    # Publication-friendly LaTeX table.
    latex_rows = sorted(
        consensus_rows,
        key=lambda r: (int(r["replication_count"]), r["signal_source"], r["signal_id"]),
        reverse=True,
    )
    columns = ["signal_source", "branch_or_site"] + pathway_names + ["replication_count", "interpretation"]
    colspec = "ll" + ("c" * len(pathway_names)) + "cl"

    latex_lines = [
        "\\begin{tabular}{" + colspec + "}",
        "\\hline",
        " & ".join(escape_latex(col) for col in columns) + r" \\",
        "\\hline",
    ]
    for row in latex_rows:
        present = set(x.strip() for x in row.get("pathways_present", "").split(",") if x.strip())
        values = [
            escape_latex(row["signal_source"]),
            escape_latex(row["signal_id"]),
        ]
        for pathway in pathway_names:
            values.append("1" if pathway in present else "0")
        values.extend([
            escape_latex(row["replication_count"]),
            escape_latex(row["interpretation"]),
        ])
        latex_lines.append(" & ".join(values) + r" \\")
    latex_lines.extend(["\\hline", "\\end{tabular}", ""])
    Path(a.latex_out).write_text("\n".join(latex_lines), encoding="utf-8")


if __name__ == "__main__":
    main()
