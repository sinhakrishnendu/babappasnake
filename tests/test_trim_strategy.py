from __future__ import annotations

from argparse import Namespace

import pytest
import yaml

from babappasnake import cli


def test_legacy_use_clipkit_mapping():
    assert cli.derive_trim_strategy(None, "yes") == "clipkit"
    assert cli.derive_trim_strategy(None, "no") == "raw"


def test_explicit_trim_strategy_overrides_legacy_flag():
    assert cli.derive_trim_strategy("raw", "yes") == "raw"
    assert cli.derive_trim_strategy("clipkit", "no") == "clipkit"
    assert cli.derive_trim_strategy("both", "no") == "both"


def test_trim_strategy_to_trim_states():
    assert cli.trim_states_from_strategy("raw") == ["raw"]
    assert cli.trim_states_from_strategy("clipkit") == ["clipkit"]
    assert cli.trim_states_from_strategy("both") == ["raw", "clipkit"]


def test_recombination_mode_alias_and_validation():
    assert cli.parse_recombination_mode("none") == "none"
    assert cli.parse_recombination_mode("gard") == "gard"
    assert cli.parse_recombination_mode("auto") == "auto"
    assert cli.effective_recombination_mode("auto") == "gard"
    assert cli.effective_recombination_mode("gard") == "gard"


def test_all_three_methods_with_both_trim_states_produces_six_unique_pathways():
    methods = cli.parse_alignment_method_option("4")
    trim_states = cli.trim_states_from_strategy("both")
    pathways = cli.enumerate_pathways(methods, trim_states)

    assert pathways == [
        "babappalign_raw",
        "babappalign_clipkit",
        "mafft_raw",
        "mafft_clipkit",
        "prank_raw",
        "prank_clipkit",
    ]
    assert len(pathways) == 6
    assert len(set(pathways)) == 6


def test_write_config_keeps_trim_strategy_and_legacy_use_clipkit_flag(tmp_path):
    outdir = tmp_path / "run"
    outdir.mkdir()

    args = Namespace(
        coverage=0.7,
        threads=12,
        orthogroup_method="rbh",
        alignment_methods="4",
        run_asr="yes",
        outgroup="culex",
        guided="yes",
        snake_args="",
        iqtree_bootstrap=1000,
        iqtree_bnni="no",
        iqtree_model="MFP",
        absrel_branches="Leaves",
        meme_branches="Leaves",
        codeml_codonfreq=7,
        clipkit_mode_protein="kpic-smart-gap",
        clipkit_mode_codon="kpic-smart-gap",
        recombination="none",
        gard_mode="Faster",
        gard_rate_classes=3,
        absrel_p=0.05,
        absrel_dynamic_start=0.05,
        absrel_dynamic_step=0.01,
        absrel_dynamic_max=0.2,
        meme_p=0.05,
    )

    cfg_path = cli.write_config(
        args=args,
        outdir=outdir,
        executables={
            "blastp": "blastp",
            "makeblastdb": "makeblastdb",
            "iqtree": "iqtree",
            "hyphy": "hyphy",
            "codeml": "codeml",
        },
        methods=["babappalign", "mafft", "prank"],
        trim_strategy="both",
        trim_states=["raw", "clipkit"],
        per_method_cores=4,
        per_pathway_cores=2,
    )

    cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))
    assert cfg["orthogroup_method"] == "rbh"
    assert cfg["trim_strategy"] == "both"
    assert cfg["trim_states"] == ["raw", "clipkit"]
    assert cfg["use_clipkit"] is True
    assert cfg["recombination"] == "none"
    assert cfg["gard_mode"] == "Faster"
    assert cfg["gard_rate_classes"] == 3
    assert cfg["per_method_cores"] == 4
    assert cfg["per_pathway_cores"] == 2
    assert cfg["alignment_method_option"] == "4"
    assert cfg["guided_mode"] is True
    assert cfg["snake_args"] == ""


def test_force_robustness_trim_strategy_always_returns_both():
    strategy, states = cli.force_robustness_trim_strategy("raw")
    assert strategy == "both"
    assert states == ["raw", "clipkit"]

    strategy, states = cli.force_robustness_trim_strategy("clipkit")
    assert strategy == "both"
    assert states == ["raw", "clipkit"]


def test_build_step_plan_includes_gard_step_when_requested():
    steps = cli.build_step_plan(
        have_cds=True,
        methods=["babappalign"],
        trim_states=["raw", "clipkit"],
        recombination_mode="gard",
        run_asr=True,
    )
    rules = [s.rule for s in steps]
    assert "gard_all_pathways" in rules

    steps_no_gard = cli.build_step_plan(
        have_cds=True,
        methods=["babappalign"],
        trim_states=["raw", "clipkit"],
        recombination_mode="none",
        run_asr=True,
    )
    rules_no_gard = [s.rule for s in steps_no_gard]
    assert "gard_all_pathways" not in rules_no_gard


def test_build_step_plan_omits_asr_steps_when_disabled():
    steps = cli.build_step_plan(
        have_cds=True,
        methods=["babappalign"],
        trim_states=["raw", "clipkit"],
        recombination_mode="none",
        run_asr=False,
    )
    rules = [s.rule for s in steps]
    assert "codeml_asr_all_pathways" not in rules
    assert "extract_selected_branch_ancestors" not in rules
    assert "asr_primary_alias" not in rules


def test_write_config_records_run_asr_flag(tmp_path):
    outdir = tmp_path / "run"
    outdir.mkdir()
    args = Namespace(
        coverage=0.7,
        threads=12,
        orthogroup_method="rbh",
        alignment_methods="1",
        run_asr="no",
        outgroup="culex",
        guided="yes",
        snake_args="",
        iqtree_bootstrap=1000,
        iqtree_bnni="no",
        iqtree_model="MFP",
        absrel_branches="Leaves",
        meme_branches="Leaves",
        codeml_codonfreq=7,
        clipkit_mode_protein="kpic-smart-gap",
        clipkit_mode_codon="kpic-smart-gap",
        recombination="none",
        gard_mode="Faster",
        gard_rate_classes=3,
        absrel_p=0.05,
        absrel_dynamic_start=0.05,
        absrel_dynamic_step=0.01,
        absrel_dynamic_max=0.2,
        meme_p=0.05,
    )
    cfg_path = cli.write_config(
        args=args,
        outdir=outdir,
        executables={
            "blastp": "blastp",
            "makeblastdb": "makeblastdb",
            "iqtree": "iqtree",
            "hyphy": "hyphy",
            "codeml": "codeml",
        },
        methods=["babappalign"],
        trim_strategy="both",
        trim_states=["raw", "clipkit"],
        per_method_cores=4,
        per_pathway_cores=2,
    )
    cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))
    assert cfg["run_asr"] is False


def test_build_snakemake_cmd_scopes_execution_to_run_directory(tmp_path):
    config_path = tmp_path / "config.yaml"
    config_path.write_text("outdir: /tmp/example\n", encoding="utf-8")
    snakefile = tmp_path / "Snakefile"
    snakefile.write_text("rule all:\n    input: []\n", encoding="utf-8")

    cmd = cli.build_snakemake_cmd(
        config_path=config_path,
        cores=4,
        target="all",
        snake_args="",
        snakefile=snakefile,
        workdir=tmp_path,
    )

    assert "--directory" in cmd
    assert str(tmp_path.resolve()) in cmd
    assert "--rerun-incomplete" in cmd


def test_prepare_resume_run_loads_saved_config_and_applies_overrides(tmp_path):
    outdir = tmp_path / "run"
    outdir.mkdir()
    config_path = outdir / "config.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "orthogroup_method": "rbh",
                "alignment_methods": ["babappalign"],
                "alignment_method_option": "1",
                "trim_strategy": "both",
                "trim_states": ["raw", "clipkit"],
                "coverage": 0.7,
                "threads": 8,
                "outgroup_query": "saved_outgroup",
                "iqtree_bootstrap": 1000,
                "iqtree_bnni": False,
                "iqtree_model": "MFP",
                "absrel_branches": "Leaves",
                "meme_branches": "Leaves",
                "codeml_codonfreq": 7,
                "run_asr": True,
                "recombination": "none",
                "gard_mode": "Faster",
                "gard_rate_classes": 3,
                "use_clipkit": True,
                "clipkit_mode_protein": "kpic-smart-gap",
                "clipkit_mode_codon": "kpic-smart-gap",
                "absrel_p": 0.05,
                "absrel_dynamic_start": 0.05,
                "absrel_dynamic_step": 0.01,
                "absrel_dynamic_max": 0.2,
                "meme_p": 0.05,
                "guided_mode": False,
                "snake_args": "--keep-going",
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    cds_path = tmp_path / "new_cds.fasta"
    cds_path.write_text(">seq\nATGAAATAG\n", encoding="utf-8")

    args = Namespace(
        prot=None,
        query=None,
        cds=str(cds_path),
        orthogroup_method="rbh",
        outdir=str(outdir),
        alignment_methods="4",
        coverage=0.7,
        threads=16,
        outgroup="override_outgroup",
        iqtree_bootstrap=1000,
        iqtree_bnni="no",
        iqtree_model="MFP",
        absrel_branches="Leaves",
        meme_branches="Leaves",
        codeml_codonfreq=7,
        run_asr="yes",
        recombination="none",
        gard_mode="Faster",
        gard_rate_classes=3,
        trim_strategy=None,
        clipkit_mode_protein="kpic-smart-gap",
        clipkit_mode_codon="kpic-smart-gap",
        absrel_p=0.05,
        absrel_dynamic_start=0.05,
        absrel_dynamic_step=0.01,
        absrel_dynamic_max=0.2,
        meme_p=0.05,
        use_clipkit="yes",
        snake_args="--latency-wait 5",
        interactive="yes",
        guided="yes",
        resume=True,
    )

    prepared, prepared_config, prepared_outdir = cli.prepare_resume_run(
        args,
        [
            "--resume",
            "--outdir",
            str(outdir),
            "--threads",
            "16",
            "--guided",
            "yes",
            "--snake-args",
            "--latency-wait 5",
            "--outgroup",
            "override_outgroup",
            "--cds",
            str(cds_path),
        ],
    )

    assert prepared_config == config_path
    assert prepared_outdir == outdir
    assert prepared.alignment_methods == "1"
    assert prepared.guided == "yes"
    assert prepared.snake_args == "--latency-wait 5"
    assert prepared.outgroup == "override_outgroup"
    assert prepared.threads == 16
    assert (outdir / "user_supplied" / "orthogroup_cds.fasta").exists()

    cfg = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    assert cfg["threads"] == 16
    assert cfg["guided_mode"] is True
    assert cfg["snake_args"] == "--latency-wait 5"
    assert cfg["outgroup_query"] == "override_outgroup"
    assert cfg["per_method_cores"] == 16
    assert cfg["per_pathway_cores"] == 8
    state = cli.load_resume_state(outdir)
    assert state["guided_mode"] is True


def test_validate_resume_request_rejects_analysis_overrides(tmp_path):
    outdir = tmp_path / "run"
    outdir.mkdir()
    (outdir / "config.yaml").write_text("orthogroup_method: rbh\n", encoding="utf-8")
    args = Namespace(outdir=str(outdir), cds=None)

    with pytest.raises(SystemExit):
        cli.validate_resume_request(args, ["--resume", "--coverage", "0.9"])


def test_validate_recombination_tools_requires_hyphy_for_gard():
    cli.validate_recombination_tools("none", {})
    cli.validate_recombination_tools("gard", {"hyphy": "/usr/bin/hyphy"})
    with pytest.raises(SystemExit):
        cli.validate_recombination_tools("gard", {})


def test_validate_orthogroup_method_tools_rbh_requires_blast_tools():
    cli.validate_orthogroup_method_tools(
        "rbh",
        {
            "blastp": "/usr/bin/blastp",
            "makeblastdb": "/usr/bin/makeblastdb",
        },
    )


def test_validate_orthogroup_method_tools_orthofinder_requires_binary():
    cli.validate_orthogroup_method_tools(
        "orthofinder",
        {
            "orthofinder": "/usr/bin/orthofinder",
            "blastp": "/usr/bin/blastp",
            "makeblastdb": "/usr/bin/makeblastdb",
        },
    )


def test_validate_orthogroup_method_tools_rbh_fallback_requires_orthofinder():
    with pytest.raises(SystemExit):
        cli.validate_orthogroup_method_tools(
            "rbh_fallback",
            {
                "blastp": "/usr/bin/blastp",
                "makeblastdb": "/usr/bin/makeblastdb",
            },
        )

    cli.validate_orthogroup_method_tools(
        "rbh_fallback",
        {
            "blastp": "/usr/bin/blastp",
            "makeblastdb": "/usr/bin/makeblastdb",
            "orthofinder": "/usr/bin/orthofinder",
        },
    )


def test_maybe_prompt_outgroup_after_cds_skips_prompt_when_cds_missing(tmp_path, monkeypatch):
    args = Namespace(outgroup="")
    config_path = tmp_path / "config.yaml"
    config_path.write_text("outgroup_query: \"\"\n", encoding="utf-8")

    prompt_calls: list[str] = []
    config_updates: list[str] = []

    def fake_prompt(*args, **kwargs):
        prompt_calls.append("called")
        return "should_not_happen"

    def fake_set(config_path_arg, outgroup_arg):
        config_updates.append(str(outgroup_arg))

    monkeypatch.setattr(cli, "prompt_text", fake_prompt)
    monkeypatch.setattr(cli, "set_outgroup_in_config", fake_set)

    cli.maybe_prompt_outgroup_after_cds(args, config_path, have_cds=False)

    assert prompt_calls == []
    assert config_updates == []
    assert args.outgroup == ""
