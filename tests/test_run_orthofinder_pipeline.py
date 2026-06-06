from __future__ import annotations

from pathlib import Path

import pytest

from babappasnake.scripts.run_orthofinder_pipeline import (
    load_member_scores_from_blast_tsv,
    load_orthogroups_from_tsv,
    parse_members,
    select_best_orthogroup_from_blast_tsv,
    select_members_for_species,
)


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def test_parse_members_handles_empty_and_comma_values():
    assert parse_members("") == []
    assert parse_members("  ") == []
    assert parse_members("a,b, c") == ["a", "b", "c"]


def test_load_orthogroups_from_tsv_returns_species_mapping(tmp_path: Path):
    table = tmp_path / "Orthogroups.tsv"
    _write(
        table,
        (
            "Orthogroup\tbabappasnake_query\tspeciesA\tspeciesB\n"
            "OG0001\tquery1\tA1, A2\tB1\n"
            "OG0002\t\tA3\t\n"
        ),
    )
    species_cols, groups = load_orthogroups_from_tsv(
        orthogroups_tsv=table,
        query_species="babappasnake_query",
    )
    assert species_cols == ["speciesA", "speciesB"]
    assert groups["OG0001"]["speciesA"] == ["A1", "A2"]
    assert groups["OG0001"]["speciesB"] == ["B1"]


def test_load_orthogroups_from_tsv_raises_if_query_species_column_missing(tmp_path: Path):
    table = tmp_path / "Orthogroups.tsv"
    _write(
        table,
        (
            "Orthogroup\tspeciesA\n"
            "OG0001\tA1\n"
        ),
    )
    with pytest.raises(RuntimeError, match="missing query species column"):
        load_orthogroups_from_tsv(
            orthogroups_tsv=table,
            query_species="babappasnake_query",
        )


def test_select_best_orthogroup_from_blast_tsv_prefers_broader_species_support(tmp_path: Path):
    blast_tsv = tmp_path / "query_vs_orthogroups.tsv"
    _write(
        blast_tsv,
        (
            "q1\tOG0002||speciesA||A2\t99.0\t90\t100\t100\t1e-40\t200\n"
            "q1\tOG0001||speciesA||A1\t98.0\t90\t100\t100\t1e-35\t180\n"
            "q1\tOG0001||speciesB||B1\t97.0\t88\t100\t100\t1e-30\t170\n"
        ),
    )
    selected = select_best_orthogroup_from_blast_tsv(blast_tsv, min_coverage=0.7)
    assert selected == "OG0001"


def test_select_best_orthogroup_from_blast_tsv_raises_on_no_passing_hits(tmp_path: Path):
    blast_tsv = tmp_path / "query_vs_orthogroups.tsv"
    _write(
        blast_tsv,
        "q1\tOG0002||speciesA||A2\t99.0\t10\t100\t100\t1e-10\t100\n",
    )
    with pytest.raises(RuntimeError, match="No BLAST hits passed coverage"):
        select_best_orthogroup_from_blast_tsv(blast_tsv, min_coverage=0.7)


def test_load_member_scores_from_blast_tsv_scores_selected_orthogroup_only(tmp_path: Path):
    blast_tsv = tmp_path / "query_vs_orthogroups.tsv"
    _write(
        blast_tsv,
        (
            "q1\tOG0001||speciesA||A1\t99.0\t90\t100\t100\t1e-35\t180\n"
            "q1\tOG0001||speciesA||A2\t98.0\t95\t100\t100\t1e-35\t220\n"
            "q1\tOG0002||speciesA||A3\t98.0\t95\t100\t100\t1e-35\t999\n"
            "q1\tOG0001||speciesA||A4\t98.0\t10\t100\t100\t1e-35\t999\n"
        ),
    )

    scores = load_member_scores_from_blast_tsv(blast_tsv, min_coverage=0.7, orthogroup_id="OG0001")

    assert scores == {"A1": 180.0, "A2": 220.0}


def test_select_members_for_species_modes_handle_multicopy_members():
    members = ["A1", "A2", "A2"]
    scores = {"A1": 100.0, "A2": 200.0}

    assert select_members_for_species(members, "strict", scores) == ([], "multi_copy_excluded")
    assert select_members_for_species(members, "representative", scores) == (["A2"], "representative_paralog")
    assert select_members_for_species(members, "paralog", scores) == (["A1", "A2"], "paralog_set")
