#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path


def read_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--query-meta", required=True)
    parser.add_argument("--threshold-summary", required=True)
    parser.add_argument("--selection-metadata", required=True)
    parser.add_argument("--members-tsv", required=True)
    parser.add_argument("--alignment-validation-json", required=True)
    parser.add_argument("--iqtree-summary-json", required=True)
    parser.add_argument("--absrel-table", required=True)
    parser.add_argument("--absrel-summary-json", required=True)
    parser.add_argument("--busted-summary-json", required=True)
    parser.add_argument("--meme-table", required=True)
    parser.add_argument("--meme-summary-json", required=True)
    parser.add_argument("--codeml-selection-json", required=True)
    parser.add_argument("--codeml-summary-tsv", required=True)
    parser.add_argument("--codeml-summary-json", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    query_meta = json.loads(Path(args.query_meta).read_text())
    thresholds = read_tsv(args.threshold_summary)
    selection = json.loads(Path(args.selection_metadata).read_text())
    members = read_tsv(args.members_tsv)
    alignment = json.loads(Path(args.alignment_validation_json).read_text())
    iqtree = json.loads(Path(args.iqtree_summary_json).read_text())
    absrel_rows = read_tsv(args.absrel_table)
    absrel_summary = json.loads(Path(args.absrel_summary_json).read_text())
    busted_summary = json.loads(Path(args.busted_summary_json).read_text())
    meme_rows = read_tsv(args.meme_table) if Path(args.meme_table).read_text().strip() else []
    meme_summary = json.loads(Path(args.meme_summary_json).read_text())
    codeml_selection = json.loads(Path(args.codeml_selection_json).read_text())
    codeml_rows = read_tsv(args.codeml_summary_tsv)
    codeml_summary = json.loads(Path(args.codeml_summary_json).read_text())

    significant_absrel = [row["branch_name"] for row in absrel_rows if row.get("significant_corrected") == "yes"]
    significant_meme = [row for row in meme_rows if row.get("significant") == "yes"][:10]
    significant_codeml = [row["branch_name"] for row in codeml_rows if row["significant_after_bh"] == "yes"]
    clipkit_enabled = bool(alignment.get("clipkit_enabled"))
    downstream_cds_is_trimmed = bool(alignment.get("downstream_cds_trimmed"))
    babappalign_protein = alignment.get("babappalign", {}).get("protein_alignment", "NA")
    babappalign_cds = alignment.get("babappalign", {}).get("cds_alignment", "NA")
    trimmed_protein = alignment.get("clipkit", {}).get("trimmed_protein_alignment", "")
    trimmed_cds = alignment.get("clipkit", {}).get("trimmed_cds_alignment", "")
    projected_trimmed_cds = alignment.get("clipkit", {}).get("projected_trimmed_cds_alignment", "")
    cds_trim_strategy = alignment.get("clipkit", {}).get("cds_trim_strategy", "untrimmed")
    effective_iqtree_alignment = alignment.get("effective_iqtree_alignment", iqtree.get("alignment", "NA"))
    effective_codon_alignment = alignment.get("effective_codon_alignment", "NA")

    lines = [
        "Executive Summary",
        "=================",
        "",
        f"Input query: {query_meta['query_id']}",
        f"Outgroup name: {selection['outgroup_name']}",
        f"Orthogroup ID: {selection['orthogroup_id']}",
        "",
        f"Selected coverage threshold: {selection['selected_threshold']}%",
        f"Selection rationale: {selection['selection_reason']}",
        "Ortholog retention by threshold:",
    ]
    for row in thresholds:
        lines.append(
            f"  - {row['threshold']}%: {row['ortholog_count']} orthologs; "
            f"outgroup_retained={row['outgroup_retained']}"
        )

    lines.extend(
        [
            "",
            f"Final orthogroup contains outgroup: {'yes' if selection['outgroup_present'] else 'no'}",
            f"Final RBH ortholog count: {selection['selected_ortholog_count']}",
            f"Final downstream analysis member count: {selection['selected_member_count']}",
            "Orthogroup members:",
        ]
    )
    for row in members:
        display_id = row.get("display_id") or row["member_id"]
        lines.append(f"  - {display_id} -> {row['member_id']} ({row['taxon_name']}, {row['member_type']})")
    if selection.get("collapsed_duplicate_count", 0):
        lines.append(
            f"Collapsed identical proteins before CDS checkpoint: {selection['collapsed_duplicate_count']}"
        )
        for row in selection.get("collapsed_duplicates", []):
            lines.append(
                f"  - kept {row['kept_member_id']} and removed {row['removed_member_id']}"
            )

    lines.extend(
        [
            "",
            f"BABAPPAlign completion: passed validation ({alignment['sequence_count']} sequences, "
            f"{alignment['alignment_length_codons']} codons in the downstream CDS alignment)",
            f"BABAPPAlign protein alignment: {babappalign_protein}",
            f"BABAPPAlign CDS alignment: {babappalign_cds}",
            f"ClipKIT enabled: {'yes' if clipkit_enabled else 'no'}",
            f"ClipKIT-trimmed protein alignment: {trimmed_protein or 'not generated'}",
            f"ClipKIT-trimmed CDS alignment: {trimmed_cds or 'not generated'}",
            f"Protein-projected CDS QC alignment: {projected_trimmed_cds or 'not generated'}",
            f"Primary CDS trim strategy: {cds_trim_strategy}",
            f"IQ-TREE mode selected by user: {iqtree.get('input_type', alignment.get('iqtree_input_type', 'NA'))}",
            f"IQ-TREE alignment used: {effective_iqtree_alignment}",
            f"IQ-TREE completion: best model {iqtree.get('best_fit_model', 'NA')}; treefile {iqtree['treefile']}",
            f"HyPhy/codeml CDS alignment used: {effective_codon_alignment}",
            "Downstream CDS alignment state: "
            f"{'trimmed' if downstream_cds_is_trimmed else 'untrimmed'}",
            "Workflow routing confirmation: IQ-TREE used the selected "
            f"{iqtree.get('input_type', alignment.get('iqtree_input_type', 'NA'))} alignment, while HyPhy/codeml "
            "used the codon-compatible CDS alignment with the same retained taxa.",
            "",
            "HyPhy findings:",
            f"  - aBSREL significant branches: {', '.join(significant_absrel) if significant_absrel else 'none'}",
            f"  - BUSTED p-value: {busted_summary.get('p_value', 'NA')} "
            f"(significant={busted_summary.get('significant', False)})",
            f"  - MEME significant sites: {meme_summary.get('significant_sites', 0)}",
        ]
    )
    if significant_meme:
        lines.append(
            "    Top MEME sites: "
            + ", ".join(f"{row['site']} (p={float(row['p_value']):.3g})" for row in significant_meme)
        )

    lines.extend(
        [
            "",
            "Selective codeml follow-up:",
            f"  - Branch selection strategy: {codeml_selection['selection_strategy_summary']}",
            "  - Tested branches:",
        ]
    )
    for branch in codeml_selection["branches"]:
        lines.append(
            f"    * {branch['branch_name']} via {branch['selection_strategy']} "
            f"(aBSREL p={branch['absrel_p_uncorrected']})"
        )

    lines.append("  - Branch-site results before/after BH correction:")
    for row in codeml_rows:
        lines.append(
            f"    * {row['branch_name']}: lnL_alt={float(row['lnL_alt']):.3f}, "
            f"lnL_null={float(row['lnL_null']):.3f}, LRT={float(row['lrt_statistic']):.3f}, "
            f"raw_p={float(row['raw_p_value']):.4g}, "
            f"bh_fdr={float(row['bh_fdr']):.4g}, significant={row['significant_after_bh']}"
        )

    interpretation = "No branch remained significant after BH correction."
    if significant_codeml:
        interpretation = (
            "Evidence for episodic positive selection remained after BH correction on: "
            + ", ".join(significant_codeml)
        )
    elif significant_absrel or busted_summary.get("significant") or meme_summary.get("significant_sites", 0) > 0:
        interpretation = (
            "Exploratory tests detected selection signals, but codeml branch-site follow-up did not retain "
            "significance after BH correction."
        )

    lines.extend(
        [
            "",
            f"Final interpretation: {interpretation}",
            "",
            "Important raw result files:",
            f"  - Selected orthogroup members: {args.members_tsv}",
            f"  - BABAPPAlign alignment validation: {args.alignment_validation_json}",
            f"  - BABAPPAlign protein alignment: {babappalign_protein}",
            f"  - BABAPPAlign CDS alignment: {babappalign_cds}",
            f"  - ClipKIT-trimmed protein alignment: {trimmed_protein or 'not generated'}",
            f"  - ClipKIT-trimmed CDS alignment: {trimmed_cds or 'not generated'}",
            f"  - Protein-projected CDS QC alignment: {projected_trimmed_cds or 'not generated'}",
            f"  - IQ-TREE summary: {args.iqtree_summary_json}",
            f"  - aBSREL parsed table: {args.absrel_table}",
            f"  - BUSTED summary: {args.busted_summary_json}",
            f"  - MEME parsed table: {args.meme_table}",
            f"  - codeml branch-site summary: {args.codeml_summary_tsv}",
            f"  - codeml branch-site JSON: {args.codeml_summary_json}",
        ]
    )

    Path(args.output).write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
