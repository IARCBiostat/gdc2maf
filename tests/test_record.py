"""The plain-text download record, and the rule that keeps it from being lost."""

import json
import time

import pandas as pd

from gdc2maf import WXS_ENSEMBLE_MAF
from gdc2maf.record import download_record_text, write_if_changed
from gdc2maf.reports import summarize_attrition


def attrition_and_files(maf_files):
    """An attrition table, a selection and a passing download check."""
    cases = pd.DataFrame(
        {
            "case_id": ["case-TCGA-AA-0001", "case-TCGA-AA-0002", "case-TCGA-AA-0009"],
            "submitter_id": ["TCGA-AA-0001", "TCGA-AA-0002", "TCGA-AA-0009"],
            "exclusion_step": pd.Series([None, None, "2. sex at birth"], dtype=object),
            "exclusion_reason": pd.Series(
                [None, None, "sex_at_birth is female"], dtype=object
            ),
        }
    )
    cases["included"] = cases["exclusion_step"].isna()
    selected = maf_files[maf_files["file_id"].isin(["f1", "f2"])].reset_index(drop=True)
    attrition = summarize_attrition("Test", cases, selected)
    check = pd.DataFrame(
        {
            "file_id": ["f1", "f2"],
            "status": ["ok", "ok"],
            "path": ["/tmp/f1", "/tmp/f2"],
        }
    )
    return attrition, selected, check


def test_record_reports_cases_excluded_and_files_downloaded(maf_files):
    attrition, selected, check = attrition_and_files(maf_files)
    text = download_record_text("Test", WXS_ENSEMBLE_MAF, attrition, selected, check)
    assert "gdc2maf download record" in text
    assert "cohort:           Test" in text
    assert "sex_at_birth is female" in text
    assert "selected (one per patient)" in text
    assert "downloaded and verified (size + md5)" in text


def test_record_names_the_gdc_release_and_client_when_recorded(tmp_path, maf_files):
    attrition, selected, check = attrition_and_files(maf_files)
    (tmp_path / "Test_query.json").write_text(
        json.dumps(
            {
                "projects": ["TCGA-LUAD", "TCGA-LUSC"],
                "sex_at_birth": ["male"],
                "gdc_data_release": "Data Release 46.0",
                "query_date": "2026-09-14",
            }
        )
    )
    (tmp_path / "Test_gdc_client.json").write_text(
        json.dumps({"version": "2.3", "source": "downloaded"})
    )
    text = download_record_text(
        "Test", WXS_ENSEMBLE_MAF, attrition, selected, check, out_dir=str(tmp_path)
    )
    assert "Data Release 46.0" in text
    assert "gdc-client:       2.3 (downloaded)" in text
    assert "TCGA-LUAD, TCGA-LUSC" in text


def test_record_holds_nothing_that_varies_between_runs(maf_files):
    attrition, selected, check = attrition_and_files(maf_files)
    first = download_record_text("Test", WXS_ENSEMBLE_MAF, attrition, selected, check)
    time.sleep(0.01)
    second = download_record_text("Test", WXS_ENSEMBLE_MAF, attrition, selected, check)
    assert first == second


def test_failed_verification_is_called_out(maf_files):
    attrition, selected, check = attrition_and_files(maf_files)
    check.loc[0, "status"] = "md5 mismatch"
    text = download_record_text("Test", WXS_ENSEMBLE_MAF, attrition, selected, check)
    assert "1  failed verification" in text
    assert "md5 mismatch: 1" in text


def test_merge_section_appears_once_the_merge_has_run(maf_files):
    attrition, selected, check = attrition_and_files(maf_files)
    sample_qc = pd.DataFrame(
        {
            "has_variants": [True, False],
            "no_variants_reason": ["", "empty GDC MAF"],
        }
    )
    text = download_record_text(
        "Test", WXS_ENSEMBLE_MAF, attrition, selected, check, sample_qc=sample_qc
    )
    assert "patients with variants" in text
    assert "empty GDC MAF: 1" in text


def test_first_write_creates_the_file(tmp_path):
    path = tmp_path / "record.txt"
    result = write_if_changed(str(path), "content\n")
    assert result["status"] == "created"
    assert path.read_text() == "content\n"


def test_rewriting_the_same_content_touches_nothing(tmp_path):
    path = tmp_path / "record.txt"
    write_if_changed(str(path), "content\n")
    before = path.stat().st_mtime_ns
    time.sleep(0.01)
    result = write_if_changed(str(path), "content\n")
    assert result["status"] == "unchanged"
    assert result["superseded"] is None
    assert path.stat().st_mtime_ns == before  # not even the mtime moved


def test_changed_content_keeps_the_previous_record(tmp_path):
    path = tmp_path / "record.txt"
    write_if_changed(str(path), "old\n")
    result = write_if_changed(str(path), "new\n")
    assert result["status"] == "updated"
    assert path.read_text() == "new\n"
    assert open(result["superseded"]).read() == "old\n"
    assert "_superseded_" in result["superseded"]


def test_several_changes_on_one_day_do_not_collide(tmp_path):
    path = tmp_path / "record.txt"
    write_if_changed(str(path), "one\n")
    first = write_if_changed(str(path), "two\n")["superseded"]
    second = write_if_changed(str(path), "three\n")["superseded"]
    assert first != second
    assert open(first).read() == "one\n"
    assert open(second).read() == "two\n"


def test_superseded_copies_can_be_turned_off(tmp_path):
    path = tmp_path / "record.txt"
    write_if_changed(str(path), "old\n")
    result = write_if_changed(str(path), "new\n", keep_superseded=False)
    assert result["superseded"] is None
    assert list(p.name for p in tmp_path.iterdir()) == ["record.txt"]


def test_logging_follows_a_stdout_replaced_after_import(tmp_path, monkeypatch):
    # configure_logging must resolve sys.stdout when called, not at import
    import io
    import logging

    from gdc2maf.provenance import configure_logging

    replacement = io.StringIO()
    monkeypatch.setattr("sys.stdout", replacement)
    configure_logging()
    try:
        logging.getLogger("gdc2maf.test").info("hello")
        assert "hello" in replacement.getvalue()
    finally:
        configure_logging(stream=None)
