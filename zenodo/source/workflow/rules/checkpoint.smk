checkpoint cds_input_checkpoint:
    input:
        members=f"{SELECTED_DIR}/orthogroup_members.tsv",
        proteins=f"{SELECTED_DIR}/orthogroup_proteins.faa",
        metadata=f"{SELECTED_DIR}/selection_metadata.json",
    output:
        directory(f"{CHECKPOINT_DIR}/request")
    log:
        f"{LOGS_DIR}/checkpoint/cds_input_checkpoint.log"
    conda:
        f"{ENV_DIR}/python.yaml"
    params:
        seed_cds=config.get("cds_input_fasta", ""),
    shell:
        r"""
        mkdir -p '{LOGS_DIR}/checkpoint'
        python '{SCRIPT_DIR}/prepare_cds_checkpoint.py' \
          --members-tsv '{input.members}' \
          --proteins-fasta '{input.proteins}' \
          --selection-metadata '{input.metadata}' \
          --checkpoint-dir '{output}' \
          --seed-cds-fasta '{params.seed_cds}' \
          > '{log}' 2>&1
        """


rule validate_user_cds:
    input:
        checkpoint_dir=cds_checkpoint_dir
    output:
        validated=f"{CHECKPOINT_DIR}/validated/cds_validated.fasta",
        member_check=f"{CHECKPOINT_DIR}/validated/cds_member_check.tsv",
        report=f"{CHECKPOINT_DIR}/validated/cds_validation_report.txt",
    log:
        f"{LOGS_DIR}/checkpoint/validate_user_cds.log"
    conda:
        f"{ENV_DIR}/python.yaml"
    shell:
        r"""
        mkdir -p '{LOGS_DIR}/checkpoint'
        python '{SCRIPT_DIR}/validate_user_cds.py' \
          --checkpoint-dir '{input.checkpoint_dir}' \
          --validated-fasta-out '{output.validated}' \
          --member-check-out '{output.member_check}' \
          --report-out '{output.report}' \
          > '{log}' 2>&1
        """
