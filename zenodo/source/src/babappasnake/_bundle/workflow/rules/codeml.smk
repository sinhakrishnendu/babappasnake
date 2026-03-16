checkpoint select_codeml_targets:
    input:
        absrel_table=f"{HYPHY_DIR}/absrel/absrel_branch_summary.tsv"
    output:
        tsv=f"{CODEML_DIR}/branch_selection.tsv",
        summary=f"{CODEML_DIR}/branch_selection.json"
    log:
        f"{LOGS_DIR}/codeml/select_codeml_targets.log"
    conda:
        f"{ENV_DIR}/python.yaml"
    params:
        fallback_count=config["codeml"]["fallback_branch_count"],
    shell:
        r"""
        mkdir -p '{CODEML_DIR}' '{LOGS_DIR}/codeml'
        python '{SCRIPT_DIR}/select_codeml_targets.py' \
          --absrel-table '{input.absrel_table}' \
          --fallback-count '{params.fallback_count}' \
          --selection-out '{output.tsv}' \
          --summary-out '{output.summary}' \
          > '{log}' 2>&1
        """


rule render_codeml_tree:
    input:
        tree=f"{IQTREE_DIR}/rooted_labeled.treefile",
        selection=lambda wildcards: checkpoints.select_codeml_targets.get().output.tsv
    output:
        treefile=f"{CODEML_DIR}/{{branch}}/foreground.treefile"
    log:
        f"{LOGS_DIR}/codeml/render_tree_{{branch}}.log"
    conda:
        f"{ENV_DIR}/python.yaml"
    shell:
        r"""
        python '{SCRIPT_DIR}/render_codeml_tree.py' \
          --treefile '{input.tree}' \
          --branch-name '{wildcards.branch}' \
          --output-tree '{output.treefile}' \
          > '{log}' 2>&1
        """


rule run_codeml_branchsite:
    input:
        alignment=effective_codon_alignment_path(),
        tree=f"{CODEML_DIR}/{{branch}}/foreground.treefile",
        selection=lambda wildcards: checkpoints.select_codeml_targets.get().output.tsv
    output:
        done=f"{CODEML_DIR}/{{branch}}/branch_site.done",
        alt_ctl=f"{CODEML_DIR}/{{branch}}/branch_site_alt/codeml.ctl",
        null_ctl=f"{CODEML_DIR}/{{branch}}/branch_site_null/codeml.ctl",
        alt_out=f"{CODEML_DIR}/{{branch}}/branch_site_alt/output.txt",
        null_out=f"{CODEML_DIR}/{{branch}}/branch_site_null/output.txt"
    log:
        f"{LOGS_DIR}/codeml/run_{{branch}}.log"
    conda:
        f"{ENV_DIR}/codeml.yaml"
    params:
        codeml=TOOL_PATHS["codeml"],
        codon_frequency=config["codeml"]["codon_frequency"],
        cleandata=config["codeml"]["cleandata"],
        kappa_initial=config["codeml"]["kappa_initial"],
        omega_initial=config["codeml"]["omega_initial"],
        noisy=config["codeml"]["noisy"],
        verbose=config["codeml"]["verbose"],
        workdir=lambda wildcards: f"{CODEML_DIR}/{wildcards.branch}",
    shell:
        r"""
        python '{SCRIPT_DIR}/run_codeml_branchsite.py' \
          --codeml '{params.codeml}' \
          --alignment-fasta '{input.alignment}' \
          --foreground-tree '{input.tree}' \
          --workdir '{params.workdir}' \
          --branch-name '{wildcards.branch}' \
          --codon-frequency '{params.codon_frequency}' \
          --cleandata '{params.cleandata}' \
          --kappa-initial '{params.kappa_initial}' \
          --omega-initial '{params.omega_initial}' \
          --noisy '{params.noisy}' \
          --verbose '{params.verbose}' \
          --done-out '{output.done}' \
          > '{log}' 2>&1
        """


rule summarize_codeml:
    input:
        selection=lambda wildcards: checkpoints.select_codeml_targets.get().output.tsv,
        done=selected_codeml_markers
    output:
        table=f"{CODEML_DIR}/codeml_branchsite_summary.tsv",
        summary=f"{CODEML_DIR}/codeml_branchsite_summary.json"
    log:
        f"{LOGS_DIR}/codeml/summarize_codeml.log"
    conda:
        f"{ENV_DIR}/codeml.yaml"
    params:
        codeml_dir=CODEML_DIR,
        fdr_alpha=config["codeml"]["fdr_alpha"],
    shell:
        r"""
        python '{SCRIPT_DIR}/summarize_codeml.py' \
          --selection-tsv '{input.selection}' \
          --codeml-dir '{params.codeml_dir}' \
          --fdr-alpha '{params.fdr_alpha}' \
          --summary-out '{output.table}' \
          --summary-json-out '{output.summary}' \
          > '{log}' 2>&1
        """
