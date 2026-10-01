"""Choosing one MAF file per patient, and saying why the others went."""

import pandas as pd
import pytest

from gdc2maf.selection import (
    SELECTION_KEYS,
    add_selection_ranks,
    select_one_file_per_case,
)


def test_ranks_prefer_primary_native_non_ffpe(maf_files):
    ranked = add_selection_ranks(maf_files).set_index("file_id")
    assert ranked.loc["f1", "rank_primary"] == 0  # 01 = primary
    assert ranked.loc["f4", "rank_primary"] == 1  # 06 = metastatic
    assert ranked.loc["f2", "rank_native"] == 0  # D = native DNA
    assert ranked.loc["f3", "rank_native"] == 1  # W = whole-genome amplified
    assert ranked.loc["f1", "rank_ffpe"] == 0
    assert ranked.loc["f1", "rank_normal"] == 0  # 10 = blood normal


def test_one_file_per_patient_is_kept(tmp_path, maf_files):
    selected, dropped = select_one_file_per_case("Test", maf_files, str(tmp_path))
    assert len(selected) == maf_files["case_id"].nunique() == 3
    assert not selected["case_id"].duplicated().any()
    assert len(dropped) == len(maf_files) - len(selected)


def test_native_dna_wins_over_wga(tmp_path, maf_files):
    selected, dropped = select_one_file_per_case("Test", maf_files, str(tmp_path))
    chosen = selected.set_index("submitter_id").loc["TCGA-AA-0002", "file_id"]
    assert chosen == "f2"
    assert dropped.set_index("file_id").loc["f3", "drop_reason"] == (
        "WGA tumour DNA (native DNA available)"
    )


def test_every_dropped_file_names_the_key_that_decided(tmp_path, maf_files):
    _, dropped = select_one_file_per_case("Test", maf_files, str(tmp_path))
    labels = {label for _, _, label in SELECTION_KEYS}
    assert set(dropped["drop_reason"]) <= labels
    assert dropped["selected_file_id"].notna().all()


def test_outputs_are_written(tmp_path, maf_files):
    select_one_file_per_case("Test", maf_files, str(tmp_path))
    assert (tmp_path / "Test_files_selected.tsv").exists()
    assert (tmp_path / "Test_files_dropped.tsv").exists()


def test_selection_is_independent_of_row_order(tmp_path, maf_files):
    forward, _ = select_one_file_per_case("A", maf_files, str(tmp_path))
    reversed_rows = maf_files.iloc[::-1].reset_index(drop=True)
    backward, _ = select_one_file_per_case("B", reversed_rows, str(tmp_path))
    assert sorted(forward["file_id"]) == sorted(backward["file_id"])


def test_identical_selection_keys_are_an_error(tmp_path, maf_files):
    twin = maf_files[maf_files["file_id"] == "f1"].copy()
    twin["file_id"] = "f1-copy"
    twin["tumor_aliquot_id"] = "tumor-f1-copy"
    doubled = pd.concat([maf_files, twin], ignore_index=True)
    with pytest.raises(ValueError, match="identical selection keys"):
        select_one_file_per_case("Test", doubled, str(tmp_path))
