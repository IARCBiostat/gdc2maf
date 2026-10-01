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
