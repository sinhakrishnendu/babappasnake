from __future__ import annotations

from argparse import Namespace

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
        outgroup="culex",
        iqtree_bootstrap=1000,
        iqtree_bnni="no",
        iqtree_model="MFP",
        absrel_branches="Leaves",
        meme_branches="Leaves",
        codeml_codonfreq=7,
        clipkit_mode_protein="kpic-smart-gap",
        clipkit_mode_codon="kpic-smart-gap",
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
    assert cfg["per_method_cores"] == 4
    assert cfg["per_pathway_cores"] == 2


def test_force_robustness_trim_strategy_always_returns_both():
    strategy, states = cli.force_robustness_trim_strategy("raw")
    assert strategy == "both"
    assert states == ["raw", "clipkit"]

    strategy, states = cli.force_robustness_trim_strategy("clipkit")
    assert strategy == "both"
    assert states == ["raw", "clipkit"]


def test_validate_orthogroup_method_tools_rbh_requires_blast_tools():
    cli.validate_orthogroup_method_tools(
        "rbh",
        {
            "blastp": "/usr/bin/blastp",
            "makeblastdb": "/usr/bin/makeblastdb",
            "orthofinder": "/usr/bin/orthofinder",
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
