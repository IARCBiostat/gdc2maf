"""The Cohort object: paths, caching, and running steps in order."""

import pandas as pd
import pytest

from gdc2maf import Cohort


def make_cohort(tmp_path, **kwargs):
    """Build a Cohort writing into ``tmp_path``."""
    return Cohort(
        name="Test",
        projects=["TCGA-TEST"],
        out_dir=str(tmp_path / "out"),
        download_dir=str(tmp_path / "downloads"),
        **kwargs,
    )


def test_out_dir_is_created(tmp_path):
    cohort = make_cohort(tmp_path)
    assert (tmp_path / "out").is_dir()


def test_output_paths_are_prefixed_with_the_cohort_name(tmp_path):
    cohort = make_cohort(tmp_path)
    assert cohort.cases_path.endswith("Test_cases.tsv")
    assert cohort.files_path.endswith("Test_files_metadata.tsv")
    assert cohort.download_check_path.endswith("Test_download_check.tsv")
    assert cohort.sample_qc_path.endswith("Test_sample_qc.tsv")
    assert cohort.log_path.endswith("Test_run.log")


def test_the_maf_is_named_after_the_file_spec(tmp_path):
    cohort = make_cohort(tmp_path)
    assert cohort.maf_path.endswith("Test_wxs_ensemble_maf.maf")


def test_an_unknown_refresh_step_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="unknown refresh step"):
        make_cohort(tmp_path, refresh=["nonsense"])


def test_known_refresh_steps_and_all_are_accepted(tmp_path):
    make_cohort(tmp_path, refresh=["cases", "download"])
    make_cohort(tmp_path, refresh=["all"])


def test_a_cached_table_is_reused_instead_of_queried(tmp_path, monkeypatch):
    cohort = make_cohort(tmp_path)
    pd.DataFrame({"case_id": ["c1"], "included": [True]}).to_csv(
        cohort.cases_path, sep="\t", index=False
    )

    def fail(*args, **kwargs):
        raise AssertionError("the GDC should not be queried when a cache exists")

    monkeypatch.setattr("gdc2maf.cohort.fetch_cases", fail)
    assert list(cohort.cases()["case_id"]) == ["c1"]


def test_refresh_bypasses_the_cache(tmp_path, monkeypatch):
    cohort = make_cohort(tmp_path, refresh=["cases"])
    pd.DataFrame({"case_id": ["stale"], "included": [True]}).to_csv(
        cohort.cases_path, sep="\t", index=False
    )
    monkeypatch.setattr(
        "gdc2maf.cohort.fetch_cases",
        lambda *a, **k: pd.DataFrame({"case_id": ["fresh"], "included": [True]}),
    )
    assert list(cohort.cases()["case_id"]) == ["fresh"]


def test_a_step_runs_only_once_per_object(tmp_path, monkeypatch):
    cohort = make_cohort(tmp_path)
    calls = []

    def once(*args, **kwargs):
        calls.append(1)
        return pd.DataFrame({"case_id": ["c1"], "included": [True]})

    monkeypatch.setattr("gdc2maf.cohort.fetch_cases", once)
    cohort.cases()
    cohort.cases()
    assert len(calls) == 1


def test_a_late_step_pulls_in_the_earlier_ones(tmp_path, monkeypatch, maf_files):
    cohort = make_cohort(tmp_path)
    monkeypatch.setattr(
        "gdc2maf.cohort.fetch_cases",
        lambda *a, **k: pd.DataFrame(
            {"case_id": maf_files["case_id"].unique(), "included": True}
        ),
    )
    monkeypatch.setattr("gdc2maf.cohort.list_maf_files", lambda *a, **k: maf_files)
    selected, dropped = cohort.selection()
    assert len(selected) == 3
    assert len(dropped) == 1


def test_quality_flags_combines_both_sources(tmp_path, monkeypatch, maf_files):
    import pandas as pd

    from gdc2maf.pancan import DO_NOT_USE_COLUMN

    cohort = make_cohort(
        tmp_path, quality_annotations_path=str(tmp_path / "quality.tsv")
    )
    (tmp_path / "quality.tsv").write_text(
        f"patient_barcode\tplatform\t{DO_NOT_USE_COLUMN}\n"
        "TCGA-AA-0001\tWXS\tTrue\n"
        "TCGA-AA-0009\tWXS\tTrue\n"  # not in the cohort
    )
    monkeypatch.setattr(
        "gdc2maf.cohort.fetch_cases",
        lambda *a, **k: pd.DataFrame({
            "case_id": maf_files["case_id"].unique(),
            "submitter_id": [
                c.replace("case-", "") for c in maf_files["case_id"].unique()
            ],
            "included": True,
        }),
    )
    monkeypatch.setattr("gdc2maf.cohort.list_maf_files", lambda *a, **k: maf_files)
    monkeypatch.setattr(
        "gdc2maf.cohort.fetch_annotations",
        lambda entities: pd.DataFrame([{
            "classification": "redaction", "category": "general",
            "case_submitter_id": "TCGA-AA-0002", "entity_id": "x",
        }]),
    )

    flags = cohort.quality_flags()
    assert set(flags["source"]) == {"GDC annotation", "PanCanAtlas"}
    # only cohort patients are flagged from the external table
    assert set(flags.loc[flags["source"] == "PanCanAtlas", "submitter_id"]) == {
        "TCGA-AA-0001"
    }
    assert cohort.quality_flags_path.endswith("Test_quality_flags.tsv")


def test_flagged_patients_are_excluded_by_default(tmp_path):
    cohort = make_cohort(tmp_path)
    assert cohort.exclude_flagged is True


def _flagging_cohort(tmp_path, monkeypatch, maf_files, **kwargs):
    """A cohort where TCGA-AA-0001 is flagged by the PanCanAtlas."""
    from gdc2maf.pancan import DO_NOT_USE_COLUMN

    (tmp_path / "quality.tsv").write_text(
        f"patient_barcode\tplatform\t{DO_NOT_USE_COLUMN}\n"
        "TCGA-AA-0001\tWXS\tTrue\n"
    )
    monkeypatch.setattr(
        "gdc2maf.cohort.fetch_cases",
        lambda *a, **k: pd.DataFrame({
            "case_id": maf_files["case_id"].unique(),
            "submitter_id": [
                c.replace("case-", "") for c in maf_files["case_id"].unique()
            ],
            "included": True,
        }),
    )
    monkeypatch.setattr("gdc2maf.cohort.list_maf_files", lambda *a, **k: maf_files)
    monkeypatch.setattr(
        "gdc2maf.cohort.fetch_annotations",
        lambda entities: pd.DataFrame(columns=["classification", "category"]),
    )
    return make_cohort(
        tmp_path, quality_annotations_path=str(tmp_path / "quality.tsv"), **kwargs
    )


def test_a_flagged_patient_is_left_out_of_the_merge(tmp_path, monkeypatch, maf_files):
    cohort = _flagging_cohort(tmp_path, monkeypatch, maf_files)
    selected, _ = cohort.selection()
    assert "TCGA-AA-0001" in set(selected["submitter_id"])

    # what merged_maf() actually hands to merge_mafs
    kept = cohort._without_flagged(selected)
    assert "TCGA-AA-0001" not in set(kept["submitter_id"])
    assert len(kept) == len(selected) - 1


def test_keeping_the_flagged_patients_merges_them(tmp_path, monkeypatch, maf_files):
    cohort = _flagging_cohort(tmp_path, monkeypatch, maf_files, exclude_flagged=False)
    selected, _ = cohort.selection()
    assert cohort._without_flagged(selected) is selected


def test_a_reused_merge_that_predates_the_exclusion_warns(
    tmp_path, monkeypatch, maf_files, caplog
):
    import logging

    cohort = _flagging_cohort(tmp_path, monkeypatch, maf_files)
    stale = pd.DataFrame({"case_id": ["case-TCGA-AA-0001"], "has_variants": [True]})
    with caplog.at_level(logging.WARNING, logger="gdc2maf.cohort"):
        cohort._warn_if_merge_predates_exclusion(stale)
    assert "refresh" in caplog.text
    assert "still contains 1" in caplog.text


def test_the_quality_table_is_downloaded_once_when_no_path_is_given(
    tmp_path, monkeypatch
):
    calls = []

    def fake_fetch(dest_dir, **kwargs):
        calls.append(dest_dir)
        return {"path": str(tmp_path / "quality.tsv")}

    monkeypatch.setattr("gdc2maf.cohort.fetch_pancan_quality_annotations", fake_fetch)
    cohort = make_cohort(tmp_path, pancan_dir=str(tmp_path / "pancan"))
    assert cohort._quality_annotations() == str(tmp_path / "quality.tsv")
    assert calls == [str(tmp_path / "pancan")]


def test_an_explicit_path_is_used_without_downloading(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "gdc2maf.cohort.fetch_pancan_quality_annotations",
        lambda *a, **k: pytest.fail("should not download"),
    )
    cohort = make_cohort(tmp_path, quality_annotations_path="given.tsv")
    assert cohort._quality_annotations() == "given.tsv"


def test_the_download_can_be_turned_off(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "gdc2maf.cohort.fetch_pancan_quality_annotations",
        lambda *a, **k: pytest.fail("should not download"),
    )
    cohort = make_cohort(tmp_path, fetch_quality_annotations=False)
    assert cohort._quality_annotations() is None


def test_a_non_tcga_cohort_does_not_download_the_pancan_table(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "gdc2maf.cohort.fetch_pancan_quality_annotations",
        lambda *a, **k: pytest.fail("the PanCanAtlas covers TCGA only"),
    )
    cohort = Cohort(
        name="Target",
        projects=["TARGET-AML"],
        out_dir=str(tmp_path / "out"),
        download_dir=str(tmp_path / "downloads"),
    )
    assert cohort._quality_annotations() is None


def test_annotations_is_a_cached_step(tmp_path):
    from gdc2maf.cohort import CACHED_STEPS

    assert "annotations" in CACHED_STEPS
    cohort = make_cohort(tmp_path, refresh=["annotations"])
    assert cohort.annotations_path.endswith("Test_gdc_annotations.tsv")
