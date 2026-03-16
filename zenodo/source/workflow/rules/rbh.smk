rule write_input_manifest:
    input:
        query=QUERY_FASTA
    output:
        manifest=f"{VALIDATION_DIR}/proteome_manifest.tsv",
        query_meta=f"{VALIDATION_DIR}/query_metadata.json",
        validation=f"{VALIDATION_DIR}/input_validation.json",
    log:
        f"{LOGS_DIR}/validation/write_input_manifest.log"
    params:
        proteome_json=PROTEOME_JSON,
        outgroup=config["outgroup_name"],
        query_taxon=config.get("query_taxon_name", ""),
        species_metadata=SPECIES_METADATA_TSV,
    shell:
        r"""
        mkdir -p '{VALIDATION_DIR}' '{LOGS_DIR}/validation'
        python '{SCRIPT_DIR}/write_input_manifest.py' \
          --query-fasta '{input.query}' \
          --outgroup-name '{params.outgroup}' \
          --query-taxon-name '{params.query_taxon}' \
          --species-metadata-tsv '{params.species_metadata}' \
          --proteome-json '{params.proteome_json}' \
          --manifest-out '{output.manifest}' \
          --query-meta-out '{output.query_meta}' \
          --validation-out '{output.validation}' \
          > '{log}' 2>&1
        """


rule run_rbh_species:
    input:
        query=QUERY_FASTA,
        manifest=f"{VALIDATION_DIR}/proteome_manifest.tsv",
        query_meta=f"{VALIDATION_DIR}/query_metadata.json",
    output:
        fwd=f"{RBH_DIR}/{{species}}/forward.tsv",
        rev=f"{RBH_DIR}/{{species}}/reverse.tsv",
        candidates=f"{RBH_DIR}/{{species}}/reciprocal_candidates.tsv",
    log:
        f"{LOGS_DIR}/rbh/{{species}}.log"
    threads: config["threads"]["rbh"]
    conda:
        f"{ENV_DIR}/rbh.yaml"
    params:
        proteome=lambda wildcards: PROTEOME_BY_ID[wildcards.species]["path"],
        taxon=lambda wildcards: PROTEOME_BY_ID[wildcards.species]["file_label"],
        search_tool=config["rbh"]["search_tool"],
        blastp=TOOL_PATHS["blastp"],
        makeblastdb=TOOL_PATHS["makeblastdb"],
        diamond=TOOL_PATHS["diamond"],
        diamond_mode=config["rbh"].get("diamond_mode", ""),
        evalue=config["rbh"]["evalue"],
        max_target_seqs=config["rbh"]["max_target_seqs"],
        workdir=lambda wildcards: f"{RBH_DIR}/{wildcards.species}",
        query_db_prefix=lambda wildcards: f"{RBH_DIR}/{wildcards.species}/query_db/query",
    shell:
        r"""
        mkdir -p '{LOGS_DIR}/rbh'
        python '{SCRIPT_DIR}/run_rbh_species.py' \
          --query-fasta '{input.query}' \
          --proteome-fasta '{params.proteome}' \
          --species-id '{wildcards.species}' \
          --taxon-name '{params.taxon}' \
          --search-tool '{params.search_tool}' \
          --blastp '{params.blastp}' \
          --makeblastdb '{params.makeblastdb}' \
          --diamond '{params.diamond}' \
          --diamond-mode '{params.diamond_mode}' \
          --evalue '{params.evalue}' \
          --threads '{threads}' \
          --max-target-seqs '{params.max_target_seqs}' \
          --query-db-prefix '{params.query_db_prefix}' \
          --workdir '{params.workdir}' \
          --forward-out '{output.fwd}' \
          --reverse-out '{output.rev}' \
          --candidate-out '{output.candidates}' \
          > '{log}' 2>&1
        """


rule compile_thresholds:
    input:
        query=QUERY_FASTA,
        query_meta=f"{VALIDATION_DIR}/query_metadata.json",
        manifest=f"{VALIDATION_DIR}/proteome_manifest.tsv",
        candidates=expand(f"{RBH_DIR}/{{species}}/reciprocal_candidates.tsv", species=SPECIES_IDS),
    output:
        summary=f"{THRESHOLD_DIR}/threshold_summary.tsv",
        summary_json=f"{THRESHOLD_DIR}/threshold_summary.json",
        members=expand(f"{THRESHOLD_DIR}/coverage_{{threshold}}/orthogroup_members.tsv", threshold=THRESHOLDS),
        fastas=expand(f"{THRESHOLD_DIR}/coverage_{{threshold}}/orthogroup_proteins.faa", threshold=THRESHOLDS),
    log:
        f"{LOGS_DIR}/rbh/compile_thresholds.log"
    conda:
        f"{ENV_DIR}/python.yaml"
    params:
        thresholds_json=json.dumps(THRESHOLDS),
        outgroup=config["outgroup_name"],
        threshold_base_dir=THRESHOLD_DIR,
    shell:
        r"""
        mkdir -p '{THRESHOLD_DIR}'
        python '{SCRIPT_DIR}/compile_thresholds.py' \
          --query-fasta '{input.query}' \
          --query-meta '{input.query_meta}' \
          --manifest '{input.manifest}' \
          --outgroup-name '{params.outgroup}' \
          --thresholds-json '{params.thresholds_json}' \
          --candidate-files {input.candidates} \
          --threshold-base-dir '{params.threshold_base_dir}' \
          --summary-out '{output.summary}' \
          --summary-json-out '{output.summary_json}' \
          > '{log}' 2>&1
        """


rule select_orthogroup:
    input:
        summary=f"{THRESHOLD_DIR}/threshold_summary.tsv",
        query_meta=f"{VALIDATION_DIR}/query_metadata.json",
        members=expand(f"{THRESHOLD_DIR}/coverage_{{threshold}}/orthogroup_members.tsv", threshold=THRESHOLDS),
        fastas=expand(f"{THRESHOLD_DIR}/coverage_{{threshold}}/orthogroup_proteins.faa", threshold=THRESHOLDS),
    output:
        metadata=f"{SELECTED_DIR}/selection_metadata.json",
        report=f"{SELECTED_DIR}/selection_report.txt",
        members=f"{SELECTED_DIR}/orthogroup_members.tsv",
        fasta=f"{SELECTED_DIR}/orthogroup_proteins.faa",
    log:
        f"{LOGS_DIR}/rbh/select_orthogroup.log"
    conda:
        f"{ENV_DIR}/python.yaml"
    params:
        threshold_base_dir=THRESHOLD_DIR,
        selected_dir=SELECTED_DIR,
        outgroup=config["outgroup_name"],
        min_ratio=config["threshold_selection"]["min_relative_retention_of_max"],
        max_loss=config["threshold_selection"]["max_absolute_loss_from_max"],
        collapse_flag="--collapse-identical-proteins"
        if config.get("post_selection", {}).get("collapse_identical_proteins", False)
        else "",
    shell:
        r"""
        mkdir -p '{SELECTED_DIR}'
        python '{SCRIPT_DIR}/select_orthogroup.py' \
          --summary-tsv '{input.summary}' \
          --threshold-base-dir '{params.threshold_base_dir}' \
          --selected-dir '{params.selected_dir}' \
          --query-meta '{input.query_meta}' \
          --outgroup-name '{params.outgroup}' \
          --min-ratio '{params.min_ratio}' \
          --max-loss '{params.max_loss}' \
          {params.collapse_flag} \
          --metadata-out '{output.metadata}' \
          --report-out '{output.report}' \
          --selected-members-out '{output.members}' \
          --selected-fasta-out '{output.fasta}' \
          > '{log}' 2>&1
        """
