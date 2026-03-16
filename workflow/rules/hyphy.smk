rule run_absrel:
    input:
        alignment=effective_codon_alignment_path(),
        tree=f"{IQTREE_DIR}/rooted_labeled.treefile",
    output:
        json=f"{HYPHY_DIR}/absrel/absrel.json",
        fit=f"{HYPHY_DIR}/absrel/absrel.fit.json",
    log:
        f"{LOGS_DIR}/hyphy/absrel.log"
    threads: config["threads"]["hyphy"]
    conda:
        f"{ENV_DIR}/hyphy.yaml"
    params:
        executable=TOOL_PATHS["hyphy"],
        branches=config["hyphy"]["branches"],
        code=config["hyphy"]["genetic_code"],
        multiple_hits=config["hyphy"]["absrel"]["multiple_hits"],
        srv=config["hyphy"]["absrel"]["srv"],
    shell:
        r"""
        mkdir -p '{HYPHY_DIR}/absrel' '{LOGS_DIR}/hyphy'
        OMP_NUM_THREADS='{threads}' '{params.executable}' absrel \
          --alignment '{input.alignment}' \
          --tree '{input.tree}' \
          --code '{params.code}' \
          --branches '{params.branches}' \
          --multiple-hits '{params.multiple_hits}' \
          --srv '{params.srv}' \
          --output '{output.json}' \
          --save-fit '{output.fit}' \
          > '{log}' 2>&1
        """


rule parse_absrel:
    input:
        json=f"{HYPHY_DIR}/absrel/absrel.json"
    output:
        table=f"{HYPHY_DIR}/absrel/absrel_branch_summary.tsv",
        summary=f"{HYPHY_DIR}/absrel/absrel_summary.json",
    log:
        f"{LOGS_DIR}/hyphy/parse_absrel.log"
    conda:
        f"{ENV_DIR}/python.yaml"
    shell:
        r"""
        python '{SCRIPT_DIR}/parse_hyphy.py' \
          --method absrel \
          --json-in '{input.json}' \
          --table-out '{output.table}' \
          --summary-out '{output.summary}' \
          > '{log}' 2>&1
        """


rule run_busted:
    input:
        alignment=effective_codon_alignment_path(),
        tree=f"{IQTREE_DIR}/rooted_labeled.treefile",
    output:
        json=f"{HYPHY_DIR}/busted/busted.json",
        fit=f"{HYPHY_DIR}/busted/busted.fit.json",
    log:
        f"{LOGS_DIR}/hyphy/busted.log"
    threads: config["threads"]["hyphy"]
    conda:
        f"{ENV_DIR}/hyphy.yaml"
    params:
        executable=TOOL_PATHS["hyphy"],
        branches=config["hyphy"]["branches"],
        code=config["hyphy"]["genetic_code"],
        multiple_hits=config["hyphy"]["busted"]["multiple_hits"],
        srv=config["hyphy"]["busted"]["srv"],
    shell:
        r"""
        mkdir -p '{HYPHY_DIR}/busted'
        OMP_NUM_THREADS='{threads}' '{params.executable}' busted \
          --alignment '{input.alignment}' \
          --tree '{input.tree}' \
          --code '{params.code}' \
          --branches '{params.branches}' \
          --multiple-hits '{params.multiple_hits}' \
          --srv '{params.srv}' \
          --output '{output.json}' \
          --save-fit '{output.fit}' \
          > '{log}' 2>&1
        """


rule parse_busted:
    input:
        json=f"{HYPHY_DIR}/busted/busted.json"
    output:
        table=f"{HYPHY_DIR}/busted/busted_table.tsv",
        summary=f"{HYPHY_DIR}/busted/busted_summary.json",
    log:
        f"{LOGS_DIR}/hyphy/parse_busted.log"
    conda:
        f"{ENV_DIR}/python.yaml"
    shell:
        r"""
        python '{SCRIPT_DIR}/parse_hyphy.py' \
          --method busted \
          --json-in '{input.json}' \
          --table-out '{output.table}' \
          --summary-out '{output.summary}' \
          > '{log}' 2>&1
        """


rule run_meme:
    input:
        alignment=effective_codon_alignment_path(),
        tree=f"{IQTREE_DIR}/rooted_labeled.treefile",
    output:
        json=f"{HYPHY_DIR}/meme/meme.json",
    log:
        f"{LOGS_DIR}/hyphy/meme.log"
    threads: config["threads"]["hyphy"]
    conda:
        f"{ENV_DIR}/hyphy.yaml"
    params:
        executable=TOOL_PATHS["hyphy"],
        branches=config["hyphy"]["branches"],
        code=config["hyphy"]["genetic_code"],
        pvalue=config["hyphy"]["meme"]["pvalue"],
        resample=config["hyphy"]["meme"]["resample"],
        multiple_hits=config["hyphy"]["meme"]["multiple_hits"],
        site_multihit=config["hyphy"]["meme"]["site_multihit"],
    shell:
        r"""
        mkdir -p '{HYPHY_DIR}/meme'
        OMP_NUM_THREADS='{threads}' '{params.executable}' meme \
          --alignment '{input.alignment}' \
          --tree '{input.tree}' \
          --code '{params.code}' \
          --branches '{params.branches}' \
          --pvalue '{params.pvalue}' \
          --resample '{params.resample}' \
          --multiple-hits '{params.multiple_hits}' \
          --site-multihit '{params.site_multihit}' \
          --output '{output.json}' \
          > '{log}' 2>&1
        """


rule parse_meme:
    input:
        json=f"{HYPHY_DIR}/meme/meme.json"
    output:
        table=f"{HYPHY_DIR}/meme/meme_sites.tsv",
        summary=f"{HYPHY_DIR}/meme/meme_summary.json",
    log:
        f"{LOGS_DIR}/hyphy/parse_meme.log"
    conda:
        f"{ENV_DIR}/python.yaml"
    params:
        threshold=config["hyphy"]["meme"]["pvalue"],
    shell:
        r"""
        python '{SCRIPT_DIR}/parse_hyphy.py' \
          --method meme \
          --json-in '{input.json}' \
          --table-out '{output.table}' \
          --summary-out '{output.summary}' \
          --meme-threshold '{params.threshold}' \
          > '{log}' 2>&1
        """
