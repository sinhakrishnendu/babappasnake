rule run_iqtree:
    input:
        alignment=effective_iqtree_alignment_path(),
        validation=f"{ALIGNMENT_DIR}/alignment_validation.json",
    output:
        treefile=f"{IQTREE_DIR}/orthogroup.treefile",
        report=f"{IQTREE_DIR}/orthogroup.iqtree",
    log:
        f"{LOGS_DIR}/iqtree/run_iqtree.log"
    threads: config["threads"]["iqtree"]
    conda:
        f"{ENV_DIR}/iqtree.yaml"
    params:
        executable=TOOL_PATHS["iqtree3"],
        prefix=f"{IQTREE_DIR}/orthogroup",
        sequence_type=effective_iqtree_sequence_type(),
        model=effective_iqtree_model(),
        bootstrap=config["iqtree"]["bootstrap"],
        alrt=config["iqtree"]["alrt"],
        extra_args=config["iqtree"].get("extra_args", ""),
        redo_flag="-redo" if config["iqtree"].get("redo", True) else "",
        input_type=IQTREE_INPUT_TYPE,
    shell:
        r"""
        mkdir -p '{IQTREE_DIR}' '{LOGS_DIR}/iqtree'
        '{params.executable}' \
          -s '{input.alignment}' \
          -pre '{params.prefix}' \
          -st '{params.sequence_type}' \
          -m '{params.model}' \
          -B '{params.bootstrap}' \
          --alrt '{params.alrt}' \
          -T '{threads}' \
          {params.redo_flag} \
          {params.extra_args} \
          > '{log}' 2>&1
        """


rule summarize_iqtree:
    input:
        report=f"{IQTREE_DIR}/orthogroup.iqtree",
        treefile=f"{IQTREE_DIR}/orthogroup.treefile",
    output:
        summary=f"{IQTREE_DIR}/iqtree_summary.json",
    log:
        f"{LOGS_DIR}/iqtree/summarize_iqtree.log"
    conda:
        f"{ENV_DIR}/python.yaml"
    params:
        alignment=effective_iqtree_alignment_path(),
        input_type=IQTREE_INPUT_TYPE,
        sequence_type=effective_iqtree_sequence_type(),
    shell:
        r"""
        python '{SCRIPT_DIR}/summarize_iqtree.py' \
          --iqtree-report '{input.report}' \
          --treefile '{input.treefile}' \
          --alignment '{params.alignment}' \
          --input-type '{params.input_type}' \
          --sequence-type '{params.sequence_type}' \
          --summary-out '{output.summary}' \
          > '{log}' 2>&1
        """


rule root_and_label_tree:
    input:
        treefile=f"{IQTREE_DIR}/orthogroup.treefile",
        members=f"{SELECTED_DIR}/orthogroup_members.tsv",
    output:
        rooted=f"{IQTREE_DIR}/rooted_labeled.treefile",
        branch_map=f"{IQTREE_DIR}/branch_map.tsv",
    log:
        f"{LOGS_DIR}/iqtree/root_and_label_tree.log"
    conda:
        f"{ENV_DIR}/python.yaml"
    params:
        outgroup=config["outgroup_name"],
    shell:
        r"""
        python '{SCRIPT_DIR}/root_and_label_tree.py' \
          --treefile '{input.treefile}' \
          --members-tsv '{input.members}' \
          --outgroup-name '{params.outgroup}' \
          --rooted-tree-out '{output.rooted}' \
          --branch-map-out '{output.branch_map}' \
          > '{log}' 2>&1
        """
