"""Attrition accounting and the one-row-per-cohort summaries."""

import pandas as pd

from gdc2maf.attrition import (
    STEP_MERGE,
    STEP_QUALITY_FLAGS,
    STEP_SEX,
    exclude_cases,
    maf_step,
)
from gdc2maf import WXS_ENSEMBLE_MAF
from gdc2maf.reports import summarize_attrition, summary_table


def cases_table():
    """Four cases: two included, one wrong sex, one without a MAF."""
    cases = pd.DataFrame(
        {
            "case_id": ["c1", "c2", "c3", "c4"],
            "submitter_id": ["p1", "p2", "p3", "p4"],
            "sex_at_birth": ["male", "male", "female", "male"],
            "has_maf": [True, True, True, False],
            "exclusion_step": pd.Series([None] * 4, dtype=object),
            "exclusion_reason": pd.Series([None] * 4, dtype=object),
        }
    )
    cases = exclude_cases(
        cases, cases["sex_at_birth"] != "male", STEP_SEX, "sex_at_birth is female"
    )
    cases = exclude_cases(
        cases, ~cases["has_maf"], maf_step(WXS_ENSEMBLE_MAF),
        f"no {WXS_ENSEMBLE_MAF.title}"
    )
    cases["included"] = cases["exclusion_step"].isna()
    return cases


def test_only_the_first_exclusion_of_a_case_is_recorded():
    cases = cases_table()
    # c3 fails both criteria but is counted once, at the sex step
    assert cases.set_index("case_id").loc["c3", "exclusion_step"] == STEP_SEX
    assert (cases["exclusion_step"].notna()).sum() == 2


def test_attrition_starts_from_every_project_case():
    cases = cases_table()
    maf_files = pd.DataFrame({"case_id": ["c1", "c2"]})
    table = summarize_attrition("Test", cases, maf_files)
    assert table.iloc[0]["n_remaining"] == 4
    assert table.iloc[-1]["n_remaining"] == 2


def test_attrition_names_the_reason_for_each_loss():
    cases = cases_table()
    maf_files = pd.DataFrame({"case_id": ["c1", "c2"]})
    reasons = set(summarize_attrition("Test", cases, maf_files)["reason"])
    assert "sex_at_birth is female" in reasons
    assert "no open WXS ensemble MAF" in reasons


def test_a_case_with_no_file_listed_is_lost_at_the_listing_step():
    cases = cases_table()
    maf_files = pd.DataFrame({"case_id": ["c1"]})  # c2 was included but has no file
    table = summarize_attrition("Test", cases, maf_files)
    assert "no MAF file returned by /files" in set(table["reason"])
    assert table.iloc[-1]["n_remaining"] == 1


def test_steps_with_no_loss_are_still_shown():
    cases = cases_table()
    maf_files = pd.DataFrame({"case_id": ["c1", "c2"]})
    table = summarize_attrition("Test", cases, maf_files)
    listing_rows = table[table["reason"] == "-"]
    assert len(listing_rows) >= 1
    assert (listing_rows["n_removed"] == 0).all()


def test_extra_steps_extend_the_table_without_the_package_knowing_them():
    cases = cases_table()
    maf_files = pd.DataFrame({"case_id": ["c1", "c2"]})
    table = summarize_attrition(
        "Test",
        cases,
        maf_files,
        extra_steps=[("7. downstream step", ["c1"], "no usable variants")],
    )
    assert "7. downstream step" in set(table["step"])
    assert "no usable variants" in set(table["reason"])
    assert table.iloc[-1]["n_remaining"] == 1


def test_sample_qc_losses_are_attributed_to_the_merge():
    cases = cases_table()
    maf_files = pd.DataFrame({"case_id": ["c1", "c2"]})
    sample_qc = pd.DataFrame(
        {
            "case_id": ["c1", "c2"],
            "has_variants": [True, False],
            "no_variants_reason": ["", "empty GDC MAF"],
        }
    )
    table = summarize_attrition("Test", cases, maf_files, sample_qc=sample_qc)
    merge_rows = table[table["step"] == STEP_MERGE]
    assert merge_rows["n_removed"].sum() == 1


def test_quality_flags_are_attributed_to_their_own_step():
    cases = cases_table()
    maf_files = pd.DataFrame({"case_id": ["c1", "c2"]})
    flags = pd.DataFrame({
        "submitter_id": ["p1", "p1"],
        "source": ["PanCanAtlas", "GDC annotation"],
        "reason": ["Do_not_use", "redaction / general"],
    })
    table = summarize_attrition("Test", cases, maf_files, quality_flags=flags)
    rows = table[table["step"] == STEP_QUALITY_FLAGS]
    assert rows["n_removed"].sum() == 1
    # both sources named, once, in one reason
    assert rows["reason"].tolist() == ["flagged by GDC annotation, PanCanAtlas"]


def test_the_quality_step_comes_before_the_merge():
    cases = cases_table()
    maf_files = pd.DataFrame({"case_id": ["c1", "c2"]})
    flags = pd.DataFrame({
        "submitter_id": ["p1"], "source": ["PanCanAtlas"], "reason": ["Do_not_use"],
    })
    sample_qc = pd.DataFrame({
        "case_id": ["c2"], "has_variants": [False],
        "no_variants_reason": ["empty GDC MAF"],
    })
    table = summarize_attrition(
        "Test", cases, maf_files, quality_flags=flags, sample_qc=sample_qc
    )
    steps = table["step"].tolist()
    assert steps.index(STEP_QUALITY_FLAGS) < steps.index(STEP_MERGE)
    # the flagged patient is gone before the merge, so the merge only sees c2
    assert table.loc[table["step"] == STEP_MERGE, "n_remaining"].iloc[-1] == 0


def test_a_patient_flagged_and_variant_less_is_counted_once():
    cases = cases_table()
    maf_files = pd.DataFrame({"case_id": ["c1", "c2"]})
    flags = pd.DataFrame({
        "submitter_id": ["p1"], "source": ["PanCanAtlas"], "reason": ["Do_not_use"],
    })
    sample_qc = pd.DataFrame({
        "case_id": ["c1"], "has_variants": [False],
        "no_variants_reason": ["empty GDC MAF"],
    })
    table = summarize_attrition(
        "Test", cases, maf_files, quality_flags=flags, sample_qc=sample_qc
    )
    assert table.loc[table["step"] == STEP_QUALITY_FLAGS, "n_removed"].sum() == 1
    assert table.loc[table["step"] == STEP_MERGE, "n_removed"].sum() == 0


def test_summary_table_stacks_cohorts_side_by_side():
    table = summary_table(
        [
            {"Cohort": "A", "Patients": 10, "Files": 12},
            {"Cohort": "B", "Patients": 20, "Files": 20},
        ]
    )
    assert list(table.columns) == ["A", "B"]
    # values are formatted as text for printing
    assert table.loc["Patients", "A"] == "10"
