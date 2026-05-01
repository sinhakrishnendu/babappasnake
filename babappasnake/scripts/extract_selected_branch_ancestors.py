#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
from dataclasses import dataclass
from io import StringIO
from pathlib import Path

from Bio import Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


NUC_CHARS = set("ACGTURYMKSWHBVDN?-")
SIGNIF_KEYS = ("significant_BH_0.05", "significant_BH_0.1", "significant_BH")


@dataclass(frozen=True)
class Pathway:
    method: str
    trim_state: str

    @property
    def name(self) -> str:
        return f"{self.method}_{self.trim_state}"


def parse_pathways(raw: str) -> list[Pathway]:
    out: list[Pathway] = []
    seen: set[tuple[str, str]] = set()
    for token in [x.strip() for x in str(raw).split(",") if x.strip()]:
        if ":" in token:
            method, trim_state = token.split(":", 1)
        elif token.endswith("_raw"):
            method, trim_state = token[:-4], "raw"
        elif token.endswith("_clipkit"):
            method, trim_state = token[:-8], "clipkit"
        else:
            method, trim_state = token, "clipkit"
        key = (method.strip(), trim_state.strip())
        if not all(key) or key in seen:
            continue
        seen.add(key)
        out.append(Pathway(*key))
    if not out:
        raise RuntimeError("No pathways were provided.")
    return out


def canonicalize(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", str(text)).strip("_")


def parse_bool(value: str) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def detect_signif_key(rows: list[dict[str, str]]) -> str | None:
    keys = {k for row in rows for k in row.keys()}
    for key in SIGNIF_KEYS:
        if key in keys:
            return key
    for key in sorted(keys):
        if key.startswith("significant_BH_"):
            return key
    return None


def hash_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def get_codeml_version(codeml_exe: str, override: str = "") -> str:
    if override.strip():
        return override.strip()
    try:
        res = subprocess.run(
            [codeml_exe],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            stdin=subprocess.DEVNULL,
        )
        text = (res.stdout or "") + "\n" + (res.stderr or "")
    except Exception as exc:
        raise RuntimeError(f"Unable to invoke codeml for version check: {exc}") from exc
    m = re.search(r"version\s+([0-9]+(?:\.[0-9A-Za-z]+)+)", text, flags=re.IGNORECASE)
    if not m:
        m = re.search(r"PAML[^0-9]*([0-9]+(?:\.[0-9A-Za-z]+)+)", text, flags=re.IGNORECASE)
    if not m:
        raise RuntimeError(
            "Could not parse codeml/PAML version from codeml output. "
            "Use --codeml-version-override if needed."
        )
    return m.group(1)


def parse_version_tuple(version_text: str) -> tuple[int, int]:
    nums = re.findall(r"\d+", version_text)
    if len(nums) < 2:
        raise RuntimeError(f"Invalid version string: {version_text}")
    return int(nums[0]), int(nums[1])


def parse_min_version(version_text: str) -> tuple[int, int]:
    m = re.match(r"^\s*(\d+)\.(\d+)\s*$", version_text)
    if not m:
        raise RuntimeError(f"Invalid --min-paml-version: {version_text}")
    return int(m.group(1)), int(m.group(2))


def read_tree(path: Path):
    return Phylo.read(str(path), "newick")


def build_parent_map(tree) -> dict[object, object | None]:
    parent: dict[object, object | None] = {tree.root: None}
    for clade in tree.find_clades(order="preorder"):
        for child in clade.clades:
            parent[child] = clade
    return parent


def assign_internal_ids(tree) -> dict[object, str]:
    ids: dict[object, str] = {}
    idx = 1
    for clade in tree.find_clades(order="preorder"):
        if not clade.is_terminal():
            ids[clade] = f"N{idx}"
            idx += 1
    return ids


def clade_signature(clade) -> tuple[str, ...]:
    leaves = []
    for tip in clade.get_terminals():
        leaves.append(canonicalize(tip.name or ""))
    return tuple(sorted(leaves))


def canonical_foreground_id(child_clade) -> str:
    if child_clade.is_terminal():
        return f"tip:{canonicalize(child_clade.name or '')}"
    sig = "|".join(clade_signature(child_clade))
    digest = hashlib.sha1(sig.encode("utf-8")).hexdigest()[:12]
    return f"lineage:{digest}"


def parse_meme_sites(meme_json: Path, threshold: float) -> set[int]:
    if not meme_json.exists() or meme_json.stat().st_size == 0:
        return set()
    payload = json.loads(meme_json.read_text(encoding="utf-8"))
    sites: set[int] = set()

    def walk(node: object) -> None:
        if isinstance(node, dict):
            site = None
            for k in ("site", "Site", "codon", "Codon", "position", "Position"):
                if k in node:
                    try:
                        site = int(str(node[k]).strip())
                    except Exception:
                        site = None
                    break
            pval = None
            for k in ("p-value", "P-value", "p", "pvalue"):
                if k in node:
                    try:
                        pval = float(node[k])
                    except Exception:
                        pval = None
                    break
            if site is not None and pval is not None and pval <= threshold:
                sites.add(site)
            for value in node.values():
                walk(value)
        elif isinstance(node, list):
            for value in node:
                walk(value)

    walk(payload)
    return sites


def parse_beb_sites(mlc_alt: Path) -> set[int]:
    if not mlc_alt.exists() or mlc_alt.stat().st_size == 0:
        return set()
    lines = mlc_alt.read_text(encoding="utf-8", errors="ignore").splitlines()
    in_beb = False
    sites: set[int] = set()
    for raw in lines:
        line = raw.strip()
        low = line.lower()
        if "bayes empirical bayes" in low:
            in_beb = True
            continue
        if not in_beb:
            continue
        if not line:
            continue
        if low.startswith("the grid") or low.startswith("naive empirical bayes"):
            continue
        if line.startswith("Positively selected sites"):
            continue
        m = re.match(r"^\s*(\d+)\s+\S+\s+([0-9]*\.?[0-9]+)(\*+)?", raw)
        if not m:
            if sites and ("time used" in low or "tree length" in low):
                break
            continue
        site = int(m.group(1))
        prob = float(m.group(2))
        star = m.group(3) or ""
        if star or prob >= 0.95:
            sites.add(site)
    return sites


def extract_numbered_tree_from_rst(rst_path: Path):
    lines = rst_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    for i, line in enumerate(lines):
        if "tree with node labels" in line.lower():
            chunks: list[str] = []
            for nxt in lines[i + 1:]:
                stripped = nxt.strip()
                if not stripped:
                    if chunks:
                        break
                    continue
                chunks.append(stripped)
                if ";" in stripped:
                    break
            candidate = "".join(chunks)
            if "(" in candidate and ";" in candidate:
                return Phylo.read(StringIO(candidate), "newick")
    return None


def parse_rst_node_sequences(rst_path: Path) -> dict[int, str]:
    lines = rst_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    in_block = False
    seqs: dict[int, str] = {}

    for raw in lines:
        line = raw.strip()
        low = line.lower()
        if "list of extant and reconstructed sequences" in low:
            in_block = True
            continue
        if not in_block:
            continue
        if not line:
            continue
        if "tree with node labels" in low and seqs:
            break

        tokens = line.split()
        seq_start = None
        for i, tok in enumerate(tokens):
            up = tok.upper()
            if up and set(up) <= NUC_CHARS and len(up) >= 3:
                seq_start = i
                break
        if seq_start is None:
            continue

        label = " ".join(tokens[:seq_start])
        seq = "".join(tokens[seq_start:]).upper()
        seq = re.sub(r"[^A-Z?-]", "", seq)
        m = re.search(r"(?:node\s*#?\s*)?(\d+)$", label, flags=re.IGNORECASE)
        if not m:
            continue
        node_id = int(m.group(1))
        seqs[node_id] = seq
    return seqs


def clade_numeric_label(clade) -> int | None:
    candidates = [getattr(clade, "name", None), getattr(clade, "confidence", None)]
    for candidate in candidates:
        if candidate is None:
            continue
        text = str(candidate).strip()
        m = re.match(r"^(\d+)(?:\.0+)?$", text)
        if m:
            return int(m.group(1))
    return None


def map_internal_ids_to_rst(tree, internal_ids: dict[object, str], rst_tree) -> dict[str, int]:
    if rst_tree is None:
        return {}

    rst_sig_to_id: dict[tuple[str, ...], int] = {}
    for clade in rst_tree.find_clades(order="preorder"):
        if clade.is_terminal():
            continue
        rst_id = clade_numeric_label(clade)
        if rst_id is None:
            continue
        sig = clade_signature(clade)
        if sig in rst_sig_to_id and rst_sig_to_id[sig] != rst_id:
            raise RuntimeError(
                f"Ambiguous rst internal-node labeling for clade signature {sig}: "
                f"{rst_sig_to_id[sig]} vs {rst_id}"
            )
        rst_sig_to_id[sig] = rst_id

    out: dict[str, int] = {}
    for clade, node_id in internal_ids.items():
        sig = clade_signature(clade)
        rst_id = rst_sig_to_id.get(sig)
        if rst_id is not None:
            out[node_id] = rst_id
    return out


def load_codon_alignment(path: Path) -> dict[str, SeqRecord]:
    recs = list(SeqIO.parse(str(path), "fasta"))
    return {rec.id: rec for rec in recs}


def resolve_pathway_tree(outdir: Path, method: str, trim_state: str) -> tuple[Path, str]:
    rooted = outdir / "tree" / method / trim_state / "orthogroup.rooted.treefile"
    if rooted.exists() and rooted.stat().st_size > 0:
        return rooted, "rooted"

    unrooted = outdir / "tree" / method / trim_state / "orthogroup.treefile"
    if unrooted.exists() and unrooted.stat().st_size > 0:
        return unrooted, "unrooted_fallback"

    return rooted, "missing"


def resolve_child_label(tree, raw_label: str):
    tip_exact = [c for c in tree.get_terminals() if (c.name or "") == raw_label]
    if len(tip_exact) == 1:
        return tip_exact[0], "tip", ""

    canon = canonicalize(raw_label)
    tip_canon = [c for c in tree.get_terminals() if canonicalize(c.name or "") == canon]
    if len(tip_canon) == 1:
        return tip_canon[0], "tip", "resolved_by_canonical_tip_name"
    if len(tip_canon) > 1:
        raise RuntimeError(
            f"Ambiguous branch label '{raw_label}' after canonical tip matching: "
            f"{', '.join(c.name or '' for c in tip_canon)}"
        )

    internal_exact = [
        c for c in tree.find_clades(order="preorder")
        if (not c.is_terminal()) and (c.name or "") == raw_label
    ]
    if len(internal_exact) == 1:
        return internal_exact[0], "internal", ""
    if len(internal_exact) > 1:
        raise RuntimeError(f"Ambiguous internal node label '{raw_label}'")

    internal_canon = [
        c for c in tree.find_clades(order="preorder")
        if (not c.is_terminal()) and canonicalize(c.name or "") == canon
    ]
    if len(internal_canon) == 1:
        return internal_canon[0], "internal", "resolved_by_canonical_internal_name"
    if len(internal_canon) > 1:
        raise RuntimeError(f"Ambiguous internal node label '{raw_label}' after canonicalization")

    raise RuntimeError(f"Could not map foreground label '{raw_label}' to a unique child node.")


def translate_codon(codon: str) -> str:
    codon = codon.upper().replace("U", "T")
    if len(codon) != 3:
        return "X"
    if "-" in codon:
        return "-"
    if any(ch not in {"A", "C", "G", "T"} for ch in codon):
        return "X"
    aa = str(Seq(codon).translate(to_stop=False)).upper()
    return aa[0] if aa else "X"


def classify_change(parent_codon: str, child_codon: str, parent_aa: str, child_aa: str) -> str:
    p = parent_codon.upper()
    c = child_codon.upper()
    if p == c:
        return "no_change"
    if len(p) != 3 or len(c) != 3 or "-" in p or "-" in c:
        return "ambiguous"
    if any(ch not in {"A", "C", "G", "T"} for ch in p + c):
        return "ambiguous"
    if parent_aa == child_aa:
        return "synonymous"
    return "nonsynonymous"


def ensure_asr(pathway: Pathway, outdir: Path, codeml: str, codonfreq: int) -> tuple[Path, Path, Path]:
    asr_dir = outdir / "asr" / pathway.method / pathway.trim_state
    rst = asr_dir / "rst"
    mlc = asr_dir / "mlc_asr.txt"
    ctl = asr_dir / "codeml_asr.ctl"
    aln = outdir / "alignments" / pathway.method / pathway.trim_state / "mapped_orthogroup_cds.analysis.fasta"
    tree, _tree_source = resolve_pathway_tree(outdir, pathway.method, pathway.trim_state)
    expected_seqfile = str(aln.resolve())
    expected_treefile = str(tree.resolve())

    def _ctl_matches() -> bool:
        if not ctl.exists():
            return False
        values: dict[str, str] = {}
        for raw in ctl.read_text(encoding="utf-8", errors="ignore").splitlines():
            if "=" not in raw:
                continue
            key, value = raw.split("=", 1)
            values[key.strip().lower()] = value.strip()
        return (
            values.get("seqfile") == expected_seqfile
            and values.get("treefile") == expected_treefile
            and str(values.get("codonfreq", "")).strip() == str(int(codonfreq))
        )

    if rst.exists() and mlc.exists() and _ctl_matches():
        return asr_dir, rst, mlc

    cmd = [
        sys.executable,
        "-m",
        "babappasnake.scripts.run_codeml_asr",
        "--alignment",
        str(aln),
        "--tree",
        str(tree),
        "--outdir",
        str(asr_dir),
        "--codeml",
        str(codeml),
        "--codonfreq",
        str(int(codonfreq)),
    ]
    res = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if res.returncode != 0 or not rst.exists() or not mlc.exists():
        raise RuntimeError(
            f"ASR auto-run failed for {pathway.name}\nCMD: {' '.join(cmd)}\n"
            f"Exit: {res.returncode}\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
        )
    return asr_dir, rst, mlc


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", required=True)
    p.add_argument("--pathways", required=True, help="Comma-separated pathway list method:trim")
    p.add_argument("--codeml", default="codeml")
    p.add_argument("--codonfreq", type=int, default=7)
    p.add_argument("--meme-p", type=float, default=0.05)
    p.add_argument("--min-paml-version", default="4.10")
    p.add_argument("--codeml-version-override", default="")
    p.add_argument("--gene", default="orthogroup")
    a = p.parse_args()

    outdir = Path(a.outdir)
    pathways = parse_pathways(a.pathways)

    min_v = parse_min_version(a.min_paml_version)
    codeml_version = get_codeml_version(str(a.codeml), str(a.codeml_version_override))
    version_tuple = parse_version_tuple(codeml_version)
    if version_tuple < min_v:
        raise RuntimeError(
            f"Detected PAML/codeml version {codeml_version} < required {a.min_paml_version}. "
            "Please use a newer PAML build."
        )

    branch_rows: list[dict[str, str]] = []
    subs_rows: list[dict[str, str]] = []
    summary_rows: list[dict[str, str]] = []
    anc_cds_records: list[SeqRecord] = []
    anc_aa_records: list[SeqRecord] = []
    desc_cds_records: list[SeqRecord] = []
    desc_aa_records: list[SeqRecord] = []

    provenance = {
        "codeml_exe": str(a.codeml),
        "codeml_version": codeml_version,
        "min_required_paml_version": str(a.min_paml_version),
        "node_labeling_scheme": "preorder_internal_ids_N1..Nn + lineage_signature_sha1",
        "pathways": [],
    }

    for pathway in pathways:
        method = pathway.method
        trim_state = pathway.trim_state
        pathway_name = pathway.name

        tree_path, tree_source = resolve_pathway_tree(outdir, method, trim_state)
        cds_aln = outdir / "alignments" / method / trim_state / "mapped_orthogroup_cds.analysis.fasta"
        branchsite_tsv = outdir / "branchsite" / method / trim_state / "branchsite_results.tsv"
        meme_json = outdir / "hyphy" / method / trim_state / "meme.json"

        asr_dir, rst_path, mlc_path = ensure_asr(pathway, outdir, str(a.codeml), int(a.codonfreq))
        ctl_path = asr_dir / "codeml_asr.ctl"

        tree = read_tree(tree_path)
        parent_map = build_parent_map(tree)
        internal_ids = assign_internal_ids(tree)
        rst_tree = extract_numbered_tree_from_rst(rst_path)
        rst_node_seqs = parse_rst_node_sequences(rst_path)
        canon_to_rst = map_internal_ids_to_rst(tree, internal_ids, rst_tree)
        codon_records = load_codon_alignment(cds_aln)
        meme_sites = parse_meme_sites(meme_json, float(a.meme_p))

        provenance["pathways"].append(
            {
                "pathway": pathway_name,
                "method": method,
                "trim_state": trim_state,
                "tree": str(tree_path),
                "tree_source": tree_source,
                "tree_sha256": hash_file(tree_path),
                "alignment": str(cds_aln),
                "alignment_sha256": hash_file(cds_aln),
                "rst": str(rst_path),
                "rst_sha256": hash_file(rst_path),
                "mlc": str(mlc_path),
                "mlc_sha256": hash_file(mlc_path),
                "ctl": str(ctl_path) if ctl_path.exists() else "",
                "ctl_sha256": hash_file(ctl_path) if ctl_path.exists() else "",
            }
        )

        bs_rows = read_tsv(branchsite_tsv)
        signif_key = detect_signif_key(bs_rows)
        selected = []
        for row in bs_rows:
            fg = str(row.get("foreground_branch", "")).strip()
            if not fg:
                continue
            if signif_key and parse_bool(str(row.get(signif_key, ""))):
                selected.append(fg)

        if not selected:
            summary_rows.append(
                {
                    "gene": a.gene,
                    "pathway": pathway_name,
                    "method": method,
                    "trim_state": trim_state,
                    "selected_branch": "NA",
                    "n_codon_changes": "0",
                    "n_nonsyn_changes": "0",
                    "n_syn_changes": "0",
                    "n_beb_overlap_changes": "0",
                    "n_meme_overlap_changes": "0",
                    "ancestral_sequence_recovered": "no",
                    "status": "no_selected_branches",
                    "notes": "",
                }
            )
            continue

        for raw_fg in selected:
            status = "ok"
            notes = ""
            try:
                safe_fg = raw_fg.replace("/", "_")
                beb_sites = parse_beb_sites(
                    outdir / "branchsite" / method / trim_state / "trees" / safe_fg / "alt" / "mlc_alt.txt"
                )
                child, child_type, resolved_note = resolve_child_label(tree, raw_fg)
                parent = parent_map.get(child)
                if parent is None:
                    raise RuntimeError("Child node has no parent (root selected as child).")
                parent_node_id = internal_ids.get(parent)
                if not parent_node_id:
                    raise RuntimeError("Parent node is not internal; cannot extract ancestral sequence.")

                if child.is_terminal():
                    child_node_id = f"TIP:{canonicalize(child.name or '')}"
                    child_name = child.name or ""
                else:
                    child_node_id = internal_ids.get(child, "NA")
                    child_name = child.name or ""

                fg_canonical = canonical_foreground_id(child)
                parent_rst_id = canon_to_rst.get(parent_node_id)
                if parent_rst_id is None:
                    raise RuntimeError(
                        f"Could not map parent internal node {parent_node_id} to an rst node id."
                    )

                parent_seq = rst_node_seqs.get(parent_rst_id, "")
                if not parent_seq:
                    raise RuntimeError(f"Missing parent ancestral sequence in rst for node {parent_rst_id}")

                if child.is_terminal():
                    child_rec = codon_records.get(child.name or "")
                    if child_rec is None:
                        # Try canonicalized lookup if exact ID fails.
                        matches = [
                            rec
                            for rec in codon_records.values()
                            if canonicalize(rec.id) == canonicalize(child.name or "")
                        ]
                        if len(matches) == 1:
                            child_rec = matches[0]
                        elif len(matches) > 1:
                            raise RuntimeError(
                                f"Ambiguous descendant sequence mapping for tip '{child.name}'."
                            )
                    if child_rec is None:
                        raise RuntimeError(f"Missing descendant codon sequence for tip '{child.name}'.")
                    child_seq = str(child_rec.seq).upper()
                else:
                    child_internal_id = internal_ids.get(child)
                    if not child_internal_id:
                        raise RuntimeError("Missing internal ID for child node.")
                    child_rst_id = canon_to_rst.get(child_internal_id)
                    if child_rst_id is None:
                        raise RuntimeError(
                            f"Could not map child internal node {child_internal_id} to an rst node id."
                        )
                    child_seq = rst_node_seqs.get(child_rst_id, "")
                    if not child_seq:
                        raise RuntimeError(
                            f"Missing child ancestral sequence in rst for node {child_rst_id}"
                        )

                parent_seq = parent_seq.upper()
                child_seq = child_seq.upper()
                n_codons = max(len(parent_seq), len(child_seq)) // 3
                n_change = 0
                n_non = 0
                n_syn = 0
                n_beb = 0
                n_meme = 0

                for i in range(n_codons):
                    start = i * 3
                    p_codon = (parent_seq[start:start + 3] if start + 3 <= len(parent_seq) else "---").upper()
                    c_codon = (child_seq[start:start + 3] if start + 3 <= len(child_seq) else "---").upper()
                    p_aa = translate_codon(p_codon)
                    c_aa = translate_codon(c_codon)
                    change_type = classify_change(p_codon, c_codon, p_aa, c_aa)
                    if change_type in {"synonymous", "nonsynonymous"}:
                        n_change += 1
                        if change_type == "synonymous":
                            n_syn += 1
                        else:
                            n_non += 1
                    overlaps_meme = "yes" if (i + 1) in meme_sites else "no"
                    overlaps_beb = "yes" if (i + 1) in beb_sites else "no"
                    if overlaps_beb == "yes" and change_type in {"synonymous", "nonsynonymous"}:
                        n_beb += 1
                    if overlaps_meme == "yes" and change_type in {"synonymous", "nonsynonymous"}:
                        n_meme += 1

                    subs_rows.append(
                        {
                            "gene": a.gene,
                            "pathway": pathway_name,
                            "method": method,
                            "trim_state": trim_state,
                            "foreground_label_canonical": fg_canonical,
                            "parent_node_id": parent_node_id,
                            "child_node_id": child_node_id,
                            "codon_site": str(i + 1),
                            "parent_codon": p_codon,
                            "child_codon": c_codon,
                            "parent_aa": p_aa,
                            "child_aa": c_aa,
                            "change_type": change_type,
                            "posterior_support_parent": "NA",
                            "posterior_support_child": "NA",
                            "overlaps_beb_site": overlaps_beb,
                            "overlaps_meme_site": overlaps_meme,
                        }
                    )

                seq_id_base = (
                    f"{a.gene}|{pathway_name}|{canonicalize(raw_fg)}|"
                    f"{parent_node_id}->{child_node_id}"
                )
                parent_aa_seq = "".join(translate_codon(parent_seq[i:i + 3]) for i in range(0, len(parent_seq) - 2, 3))
                child_aa_seq = "".join(translate_codon(child_seq[i:i + 3]) for i in range(0, len(child_seq) - 2, 3))
                anc_cds_records.append(SeqRecord(Seq(parent_seq), id=seq_id_base + "|ancestor_cds", description=""))
                desc_cds_records.append(SeqRecord(Seq(child_seq), id=seq_id_base + "|descendant_cds", description=""))
                anc_aa_records.append(SeqRecord(Seq(parent_aa_seq), id=seq_id_base + "|ancestor_aa", description=""))
                desc_aa_records.append(SeqRecord(Seq(child_aa_seq), id=seq_id_base + "|descendant_aa", description=""))

                branch_rows.append(
                    {
                        "gene": a.gene,
                        "pathway": pathway_name,
                        "method": method,
                        "trim_state": trim_state,
                        "foreground_label_raw": raw_fg,
                        "foreground_label_canonical": fg_canonical,
                        "parent_node_id": parent_node_id,
                        "child_node_id": child_node_id,
                        "child_node_name": child_name,
                        "child_node_type": child_type,
                        "status": "ok",
                        "notes": resolved_note,
                    }
                )
                summary_rows.append(
                    {
                        "gene": a.gene,
                        "pathway": pathway_name,
                        "method": method,
                        "trim_state": trim_state,
                        "selected_branch": raw_fg,
                        "n_codon_changes": str(n_change),
                        "n_nonsyn_changes": str(n_non),
                        "n_syn_changes": str(n_syn),
                        "n_beb_overlap_changes": str(n_beb),
                        "n_meme_overlap_changes": str(n_meme),
                        "ancestral_sequence_recovered": "yes",
                        "status": "ok",
                        "notes": resolved_note,
                    }
                )
            except Exception as exc:
                status = "failed"
                notes = str(exc)
                branch_rows.append(
                    {
                        "gene": a.gene,
                        "pathway": pathway_name,
                        "method": method,
                        "trim_state": trim_state,
                        "foreground_label_raw": raw_fg,
                        "foreground_label_canonical": canonicalize(raw_fg),
                        "parent_node_id": "NA",
                        "child_node_id": "NA",
                        "child_node_name": "",
                        "child_node_type": "unknown",
                        "status": status,
                        "notes": notes,
                    }
                )
                summary_rows.append(
                    {
                        "gene": a.gene,
                        "pathway": pathway_name,
                        "method": method,
                        "trim_state": trim_state,
                        "selected_branch": raw_fg,
                        "n_codon_changes": "0",
                        "n_nonsyn_changes": "0",
                        "n_syn_changes": "0",
                        "n_beb_overlap_changes": "0",
                        "n_meme_overlap_changes": "0",
                        "ancestral_sequence_recovered": "no",
                        "status": status,
                        "notes": notes,
                    }
                )

    asr_dir = outdir / "asr"
    asr_dir.mkdir(parents=True, exist_ok=True)

    write_tsv(
        asr_dir / "branch_to_nodes.tsv",
        [
            "gene",
            "pathway",
            "method",
            "trim_state",
            "foreground_label_raw",
            "foreground_label_canonical",
            "parent_node_id",
            "child_node_id",
            "child_node_name",
            "child_node_type",
            "status",
            "notes",
        ],
        branch_rows,
    )
    write_tsv(
        asr_dir / "branch_substitutions.tsv",
        [
            "gene",
            "pathway",
            "method",
            "trim_state",
            "foreground_label_canonical",
            "parent_node_id",
            "child_node_id",
            "codon_site",
            "parent_codon",
            "child_codon",
            "parent_aa",
            "child_aa",
            "change_type",
            "posterior_support_parent",
            "posterior_support_child",
            "overlaps_beb_site",
            "overlaps_meme_site",
        ],
        subs_rows,
    )
    write_tsv(
        asr_dir / "selected_branch_asr_summary.tsv",
        [
            "gene",
            "pathway",
            "method",
            "trim_state",
            "selected_branch",
            "n_codon_changes",
            "n_nonsyn_changes",
            "n_syn_changes",
            "n_beb_overlap_changes",
            "n_meme_overlap_changes",
            "ancestral_sequence_recovered",
            "status",
            "notes",
        ],
        summary_rows,
    )

    SeqIO.write(anc_cds_records, asr_dir / "ancestor_sequences_cds.fasta", "fasta")
    SeqIO.write(anc_aa_records, asr_dir / "ancestor_sequences_aa.fasta", "fasta")
    SeqIO.write(desc_cds_records, asr_dir / "descendant_sequences_cds.fasta", "fasta")
    SeqIO.write(desc_aa_records, asr_dir / "descendant_sequences_aa.fasta", "fasta")
    (asr_dir / "asr_extraction_provenance.json").write_text(
        json.dumps(provenance, indent=2),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
