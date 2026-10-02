"""Resolving a cohort into cases, and what the sex criterion does."""

import pandas as pd
import pytest

from gdc2maf import WXS_ENSEMBLE_MAF
from gdc2maf.attrition import STEP_SEX, maf_step
from gdc2maf.cases import fetch_cases, fill_sex_from_table


@pytest.fixture
def no_network(monkeypatch):
    """Answer the three GDC calls fetch_cases makes with fixed data."""
    project_cases = pd.DataFrame(
        {
            "case_id": ["c1", "c2", "c3", "c4"],
            "submitter_id": ["TCGA-AA-0001", "TCGA-AA-0002", "TCGA-AA-0003",
                             "TCGA-AA-0004"],
            "project_id": ["TCGA-TEST"] * 4,
            "primary_site": ["lung"] * 4,
            "disease_type": ["adenomas"] * 4,
            "sex_at_birth": ["male", "female", None, "male"],
        }
    )
    monkeypatch.setattr("gdc2maf.cases.gdc_cohort_cases", lambda **k: project_cases)
    # c4 has no matching MAF
    monkeypatch.setattr(
        "gdc2maf.cases.cases_with_files", lambda spec, ids: {"c1", "c2", "c3"}
    )
    monkeypatch.setattr("gdc2maf.cases.gdc_data_release", lambda: "Data Release 46.0")
    return project_cases


def test_no_sex_criterion_keeps_every_sex_including_unknown(no_network, tmp_path):
    cases = fetch_cases("Test", ["TCGA-TEST"], str(tmp_path), sex_at_birth=None)
    assert STEP_SEX not in set(cases["exclusion_step"].dropna())
    # all three with a MAF are kept, whatever their sex (c3's is missing)
    assert set(cases.loc[cases["included"], "case_id"]) == {"c1", "c2", "c3"}


def test_a_sex_criterion_excludes_the_others_and_the_unknowns(no_network, tmp_path):
    cases = fetch_cases("Test", ["TCGA-TEST"], str(tmp_path), sex_at_birth=["male"])
    included = set(cases.loc[cases["included"], "case_id"])
    assert included == {"c1"}
    reasons = dict(zip(cases["case_id"], cases["exclusion_reason"]))
    assert reasons["c2"] == "sex_at_birth is female"
    assert reasons["c3"] == "sex_at_birth missing in GDC"


def test_cases_without_a_matching_file_are_excluded_by_the_spec(no_network, tmp_path):
    cases = fetch_cases("Test", ["TCGA-TEST"], str(tmp_path), sex_at_birth=None)
    row = cases.set_index("case_id").loc["c4"]
    assert row["exclusion_step"] == maf_step(WXS_ENSEMBLE_MAF)
    assert row["exclusion_reason"] == f"no {WXS_ENSEMBLE_MAF.title}"


def test_every_project_case_is_kept_in_the_table(no_network, tmp_path):
    cases = fetch_cases("Test", ["TCGA-TEST"], str(tmp_path), sex_at_birth=["male"])
    assert len(cases) == 4  # nothing is dropped, only marked


def test_the_query_record_is_written(no_network, tmp_path):
    import json

    fetch_cases("Test", ["TCGA-TEST"], str(tmp_path), sex_at_birth=None)
    record = json.loads((tmp_path / "Test_query.json").read_text())
    assert record["projects"] == ["TCGA-TEST"]
    assert record["sex_at_birth"] is None
    assert record["gdc_data_release"] == "Data Release 46.0"
    assert record["file_spec"]["name"] == "wxs_ensemble_maf"


def test_sex_fallback_fills_only_what_the_gdc_lacks(tmp_path):
    cases = pd.DataFrame(
        {
            "submitter_id": ["TCGA-AA-0001", "TCGA-AA-0003"],
            "sex_at_birth": ["male", None],
        }
    )
    table = tmp_path / "pancan.tsv"
    table.write_text(
        "bcr_patient_barcode\tgender\nTCGA-AA-0001\tFEMALE\nTCGA-AA-0003\tFEMALE\n"
    )
    filled = fill_sex_from_table(cases, str(table))
    # the GDC value wins where both have one
    assert list(filled["sex_at_birth"]) == ["male", "female"]
    assert list(filled["sex_at_birth_source"]) == ["GDC", "PanCan"]
