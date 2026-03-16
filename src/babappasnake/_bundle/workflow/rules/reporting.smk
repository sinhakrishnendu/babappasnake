rule executive_summary:
    input:
        query_meta=f"{VALIDATION_DIR}/query_metadata.json",
        threshold_summary=f"{THRESHOLD_DIR}/threshold_summary.tsv",
        selection_metadata=f"{SELECTED_DIR}/selection_metadata.json",
        members=f"{SELECTED_DIR}/orthogroup_members.tsv",
        alignment_validation=f"{ALIGNMENT_DIR}/alignment_validation.json",
        iqtree_summary=f"{IQTREE_DIR}/iqtree_summary.json",
        absrel_table=f"{HYPHY_DIR}/absrel/absrel_branch_summary.tsv",
        absrel_summary=f"{HYPHY_DIR}/absrel/absrel_summary.json",
        busted_summary=f"{HYPHY_DIR}/busted/busted_summary.json",
        meme_table=f"{HYPHY_DIR}/meme/meme_sites.tsv",
        meme_summary=f"{HYPHY_DIR}/meme/meme_summary.json",
        codeml_selection=f"{CODEML_DIR}/branch_selection.json",
        codeml_table=f"{CODEML_DIR}/codeml_branchsite_summary.tsv",
        codeml_summary=f"{CODEML_DIR}/codeml_branchsite_summary.json",
    output:
        f"{REPORT_DIR}/executive_summary.txt"
    log:
        f"{LOGS_DIR}/reports/executive_summary.log"
    conda:
        f"{ENV_DIR}/python.yaml"
    shell:
        r"""
        mkdir -p '{REPORT_DIR}' '{LOGS_DIR}/reports'
        python '{SCRIPT_DIR}/write_executive_summary.py' \
          --query-meta '{input.query_meta}' \
          --threshold-summary '{input.threshold_summary}' \
          --selection-metadata '{input.selection_metadata}' \
          --members-tsv '{input.members}' \
          --alignment-validation-json '{input.alignment_validation}' \
          --iqtree-summary-json '{input.iqtree_summary}' \
          --absrel-table '{input.absrel_table}' \
          --absrel-summary-json '{input.absrel_summary}' \
          --busted-summary-json '{input.busted_summary}' \
          --meme-table '{input.meme_table}' \
          --meme-summary-json '{input.meme_summary}' \
          --codeml-selection-json '{input.codeml_selection}' \
          --codeml-summary-tsv '{input.codeml_table}' \
          --codeml-summary-json '{input.codeml_summary}' \
          --output '{output}' \
          > '{log}' 2>&1
        """
