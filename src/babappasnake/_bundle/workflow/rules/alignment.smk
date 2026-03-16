rule run_babappalign:
    input:
        cds=f"{CHECKPOINT_DIR}/validated/cds_validated.fasta"
    output:
        staged_input=f"{ALIGNMENT_BAB_DIR}/orthogroup_cds_input.fasta",
        raw_proteins=f"{ALIGNMENT_BAB_DIR}/orthogroup_cds_input.protein.aln.fasta",
        raw_cds=f"{ALIGNMENT_BAB_DIR}/orthogroup_cds_input.codon.aln.fasta",
        aligned_proteins=f"{ALIGNMENT_BAB_PROTEIN_DIR}/aligned_proteins.fasta",
        aligned_cds=f"{ALIGNMENT_BAB_CDS_DIR}/aligned_cds.fasta",
    log:
        f"{LOGS_DIR}/alignment/babappalign.log"
    conda:
        f"{ENV_DIR}/babappalign.yaml"
    params:
        executable=TOOL_PATHS["babappalign"],
        model_path=BABAPPALIGN_MODEL,
        device=config["alignment"]["device"],
        gap_open=config["alignment"]["gap_open"],
        gap_extend=config["alignment"]["gap_extend"],
        output_dir=ALIGNMENT_BAB_DIR,
    shell:
        r"""
        mkdir -p '{ALIGNMENT_BAB_DIR}' '{ALIGNMENT_BAB_PROTEIN_DIR}' '{ALIGNMENT_BAB_CDS_DIR}' '{LOGS_DIR}/alignment'
        python '{SCRIPT_DIR}/run_babappalign.py' \
          --executable '{params.executable}' \
          --model-path '{params.model_path}' \
          --device '{params.device}' \
          --gap-open '{params.gap_open}' \
          --gap-extend '{params.gap_extend}' \
          --input-fasta '{input.cds}' \
          --output-dir '{params.output_dir}' \
          --staged-input-out '{output.staged_input}' \
          --raw-cds-out '{output.raw_cds}' \
          --raw-proteins-out '{output.raw_proteins}' \
          --aligned-cds-out '{output.aligned_cds}' \
          --aligned-proteins-out '{output.aligned_proteins}' \
          > '{log}' 2>&1
        """


rule run_clipkit_protein:
    input:
        alignment=f"{ALIGNMENT_BAB_PROTEIN_DIR}/aligned_proteins.fasta"
    output:
        trimmed=f"{ALIGNMENT_CLIPKIT_PROTEIN_DIR}/trimmed_proteins.fasta",
        clipkit_log=f"{ALIGNMENT_CLIPKIT_PROTEIN_DIR}/trimmed_proteins.clipkit.log",
    log:
        f"{LOGS_DIR}/alignment/clipkit_protein.log"
    conda:
        f"{ENV_DIR}/clipkit.yaml"
    params:
        executable=TOOL_PATHS["clipkit"],
        mode=config["clipkit"]["mode"],
        seq_type="aa",
        threads=config["threads"].get("clipkit", 1),
    shell:
        r"""
        mkdir -p '{ALIGNMENT_CLIPKIT_PROTEIN_DIR}' '{LOGS_DIR}/alignment'
        python '{SCRIPT_DIR}/run_clipkit.py' \
          --executable '{params.executable}' \
          --input-fasta '{input.alignment}' \
          --output-fasta '{output.trimmed}' \
          --clipkit-log-out '{output.clipkit_log}' \
          --mode '{params.mode}' \
          --sequence-type '{params.seq_type}' \
          --threads '{params.threads}' \
          > '{log}' 2>&1
        """


rule run_clipkit_cds_direct:
    input:
        alignment=f"{ALIGNMENT_BAB_CDS_DIR}/aligned_cds.fasta"
    output:
        trimmed=f"{ALIGNMENT_CLIPKIT_CDS_DIR}/trimmed_cds.fasta",
        clipkit_log=f"{ALIGNMENT_CLIPKIT_CDS_DIR}/trimmed_cds.clipkit.log",
    log:
        f"{LOGS_DIR}/alignment/clipkit_cds.log"
    conda:
        f"{ENV_DIR}/clipkit.yaml"
    params:
        executable=TOOL_PATHS["clipkit"],
        mode=config["clipkit"]["mode"],
        seq_type="nt",
        threads=config["threads"].get("clipkit", 1),
    shell:
        r"""
        mkdir -p '{ALIGNMENT_CLIPKIT_CDS_DIR}' '{LOGS_DIR}/alignment'
        python '{SCRIPT_DIR}/run_clipkit.py' \
          --executable '{params.executable}' \
          --input-fasta '{input.alignment}' \
          --output-fasta '{output.trimmed}' \
          --clipkit-log-out '{output.clipkit_log}' \
          --mode '{params.mode}' \
          --sequence-type '{params.seq_type}' \
          --threads '{params.threads}' \
          --codon \
          > '{log}' 2>&1
        """


rule project_clipkit_mask_to_cds:
    input:
        protein_alignment=f"{ALIGNMENT_BAB_PROTEIN_DIR}/aligned_proteins.fasta",
        cds_alignment=f"{ALIGNMENT_BAB_CDS_DIR}/aligned_cds.fasta",
        trimmed_protein=f"{ALIGNMENT_CLIPKIT_PROTEIN_DIR}/trimmed_proteins.fasta",
        protein_clipkit_log=f"{ALIGNMENT_CLIPKIT_PROTEIN_DIR}/trimmed_proteins.clipkit.log",
        direct_cds=f"{ALIGNMENT_CLIPKIT_CDS_DIR}/trimmed_cds.fasta",
    output:
        trimmed_cds=f"{ALIGNMENT_CLIPKIT_CDS_DIR}/projected_trimmed_cds.fasta",
        report=f"{ALIGNMENT_CLIPKIT_CDS_DIR}/projection_report.txt",
        summary=f"{ALIGNMENT_CLIPKIT_CDS_DIR}/projection_summary.json",
    log:
        f"{LOGS_DIR}/alignment/project_clipkit_mask_to_cds.log"
    conda:
        f"{ENV_DIR}/python.yaml"
    shell:
        r"""
        python '{SCRIPT_DIR}/project_protein_trim_to_cds.py' \
          --protein-alignment '{input.protein_alignment}' \
          --cds-alignment '{input.cds_alignment}' \
          --trimmed-protein-alignment '{input.trimmed_protein}' \
          --protein-clipkit-log '{input.protein_clipkit_log}' \
          --direct-cds-clipkit-alignment '{input.direct_cds}' \
          --trimmed-cds-out '{output.trimmed_cds}' \
          --report-out '{output.report}' \
          --summary-out '{output.summary}' \
          > '{log}' 2>&1
        """


rule validate_alignment:
    input:
        alignment_validation_dependencies
    output:
        json=f"{ALIGNMENT_DIR}/alignment_validation.json",
        report=f"{ALIGNMENT_DIR}/alignment_validation.txt",
    log:
        f"{LOGS_DIR}/alignment/validate_alignment.log"
    conda:
        f"{ENV_DIR}/python.yaml"
    params:
        babappalign_proteins=babappalign_protein_alignment(),
        babappalign_cds=babappalign_cds_alignment(),
        trimmed_proteins=clipkit_trimmed_protein_alignment() if CLIPKIT_ENABLED else "",
        trimmed_cds=clipkit_trimmed_cds_alignment() if CLIPKIT_ENABLED else "",
        projected_trimmed_cds=projected_clipkit_trimmed_cds_alignment() if CLIPKIT_ENABLED else "",
        projection_summary=f"{ALIGNMENT_CLIPKIT_CDS_DIR}/projection_summary.json" if CLIPKIT_ENABLED else "",
        members=f"{SELECTED_DIR}/orthogroup_members.tsv",
        clipkit_enabled=str(CLIPKIT_ENABLED).lower(),
        iqtree_input_type=IQTREE_INPUT_TYPE,
        effective_iqtree_alignment=effective_iqtree_alignment_path(),
        effective_codon_alignment=effective_codon_alignment_path(),
    shell:
        r"""
        python '{SCRIPT_DIR}/validate_alignment_products.py' \
          --babappalign-protein '{params.babappalign_proteins}' \
          --babappalign-cds '{params.babappalign_cds}' \
          --trimmed-protein '{params.trimmed_proteins}' \
          --trimmed-cds '{params.trimmed_cds}' \
          --projected-trimmed-cds '{params.projected_trimmed_cds}' \
          --projection-summary '{params.projection_summary}' \
          --members-tsv '{params.members}' \
          --clipkit-enabled '{params.clipkit_enabled}' \
          --iqtree-input-type '{params.iqtree_input_type}' \
          --effective-iqtree-alignment '{params.effective_iqtree_alignment}' \
          --effective-codon-alignment '{params.effective_codon_alignment}' \
          --json-out '{output.json}' \
          --report-out '{output.report}' \
          > '{log}' 2>&1
        """
