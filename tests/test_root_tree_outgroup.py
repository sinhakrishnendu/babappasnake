from __future__ import annotations

from pathlib import Path

from babappasnake.scripts.root_tree_outgroup import apply_outgroup_or_copy


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def test_no_outgroup_copies_tree(tmp_path):
    tree = tmp_path / "orthogroup.treefile"
    output = tmp_path / "orthogroup.rooted.treefile"
    _write(tree, "(A:0.1,B:0.2,C:0.3);\n")

    status = apply_outgroup_or_copy(tree, output, "")

    assert status["status"] == "copied_unrooted"
    assert output.read_text(encoding="utf-8") == tree.read_text(encoding="utf-8")


def test_unmatched_outgroup_falls_back_to_unrooted_tree(tmp_path, capsys):
    tree = tmp_path / "orthogroup.treefile"
    output = tmp_path / "orthogroup.rooted.treefile"
    _write(tree, "(A:0.1,B:0.2,C:0.3);\n")

    status = apply_outgroup_or_copy(tree, output, "culex")
    captured = capsys.readouterr()

    assert status["status"] == "copied_unrooted"
    assert "using unrooted tree downstream" in captured.err
    assert output.read_text(encoding="utf-8") == tree.read_text(encoding="utf-8")


def test_overbroad_outgroup_falls_back_to_unrooted_tree(tmp_path, capsys):
    tree = tmp_path / "orthogroup.treefile"
    output = tmp_path / "orthogroup.rooted.treefile"
    _write(tree, "(speciesA:0.1,speciesB:0.2,speciesC:0.3);\n")

    status = apply_outgroup_or_copy(tree, output, "species")
    captured = capsys.readouterr()

    assert status["status"] == "copied_unrooted"
    assert "using unrooted tree downstream" in captured.err
    assert output.read_text(encoding="utf-8") == tree.read_text(encoding="utf-8")
