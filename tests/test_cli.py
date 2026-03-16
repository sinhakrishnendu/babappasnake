from pathlib import Path

import pytest
import yaml

from babappasnake.cli import build_parser, create_runtime_config


def test_cli_rejects_invalid_clipkit_value():
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "--input",
                "proteomes",
                "--query",
                "query.fasta",
                "--outgroup",
                "Culex quinquefasciatus",
                "--clipkit",
                "maybe",
                "--iqtree",
                "protein",
            ]
        )


def test_cli_rejects_invalid_iqtree_value():
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "--input",
                "proteomes",
                "--query",
                "query.fasta",
                "--outgroup",
                "Culex quinquefasciatus",
                "--clipkit",
                "yes",
                "--iqtree",
                "peptide",
            ]
        )


def test_runtime_config_writes_clipkit_and_iqtree_settings(tmp_path):
    parser = build_parser()
    proteomes = tmp_path / "proteomes"
    proteomes.mkdir()
    query = tmp_path / "query.fasta"
    query.write_text(">q\nMKT\n")
    base_config = tmp_path / "base.yaml"
    base_config.write_text("results_dir: default_results\nclipkit:\n  enabled: true\niqtree:\n  input_type: protein\n")

    args = parser.parse_args(
        [
            "--input",
            str(proteomes),
            "--query",
            str(query),
            "--outgroup",
            "Culex quinquefasciatus",
            "--clipkit",
            "no",
            "--iqtree",
            "cds",
            "--output",
            str(tmp_path / "results"),
            "--config",
            str(base_config),
        ]
    )
    config, runtime_path = create_runtime_config(args, Path("/unused"))

    assert runtime_path.exists()
    written = yaml.safe_load(runtime_path.read_text())
    assert config["clipkit"]["enabled"] is False
    assert written["clipkit"]["enabled"] is False
    assert written["iqtree"]["input_type"] == "cds"
