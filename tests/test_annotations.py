"""GDC annotations, and which of them question the material."""

import pandas as pd

from gdc2maf.annotations import (
    EXCLUSION_CATEGORIES,
    EXCLUSION_CLASSIFICATIONS,
    exclusion_flags,
    fetch_annotations,
    summarize_annotations,
)


def annotation(classification, category, case="TCGA-AA-0001", entity="e1"):
    """One annotation row as fetch_annotations returns it."""
    return {
        "annotation_id": f"a-{entity}-{category}",
        "entity_type": "aliquot",
        "entity_id": entity,
        "entity_submitter_id": f"{case}-01A-11D-A271-08",
        "classification": classification,
        "category": category,
        "status": "Approved",
        "notes": "",
        "created_datetime": "2024-01-01",
        "case_id": f"case-{case}",
        "case_submitter_id": case,
        "project_id": "TCGA-TEST",
    }


def test_no_entities_needs_no_query(monkeypatch):
    def fail(*a, **k):
        raise AssertionError("the GDC should not be queried for an empty entity list")

    monkeypatch.setattr("gdc2maf.annotations.gdc_query", fail)
    assert fetch_annotations([]).empty


def test_entities_are_deduplicated_and_nulls_dropped(monkeypatch):
    seen = {}

    def capture(endpoint, filters, fields):
        seen["values"] = filters["content"]["value"]
        return []

    monkeypatch.setattr("gdc2maf.annotations.gdc_query", capture)
    fetch_annotations(["a", "a", None, "b", float("nan")])
    assert seen["values"] == ["a", "b"]


def test_the_somatic_sniper_notification_is_not_exclusion_grade():
    # every open WXS ensemble MAF carries this one; treating an annotation's
    # presence as a verdict would drop the whole cohort
    found = pd.DataFrame([annotation("notification", "general")])
    assert exclusion_flags(found).empty


def test_clinical_history_notifications_are_not_exclusion_grade():
    kept = ["prior malignancy", "neoadjuvant therapy", "synchronous malignancy",
            "item is noncanonical", "bcr notification"]
    found = pd.DataFrame([annotation("notification", c) for c in kept])
    assert exclusion_flags(found).empty


def test_a_redaction_is_exclusion_grade():
    found = pd.DataFrame([annotation("redaction", "general")])
    assert len(exclusion_flags(found)) == 1


def test_the_do_not_use_category_is_exclusion_grade():
    found = pd.DataFrame([annotation("centernotification", "item flagged dnu")])
    assert len(exclusion_flags(found)) == 1


def test_identity_and_integrity_categories_are_exclusion_grade():
    for category in EXCLUSION_CATEGORIES:
        found = pd.DataFrame([annotation("notification", category)])
        assert len(exclusion_flags(found)) == 1, category


def test_the_exclusion_sets_are_deliberately_narrow():
    assert EXCLUSION_CLASSIFICATIONS == ("redaction",)
    assert "general" not in EXCLUSION_CATEGORIES
    assert "prior malignancy" not in EXCLUSION_CATEGORIES


def test_the_sets_can_be_widened_by_the_caller():
    found = pd.DataFrame([annotation("notification", "prior malignancy")])
    assert len(exclusion_flags(found, categories=["prior malignancy"])) == 1


def test_empty_annotations_summarise_without_error():
    empty = pd.DataFrame()
    summary = summarize_annotations("Test", empty, empty)
    assert summary["Annotations found"] == 0
    assert summary["Patients with an exclusion-grade annotation"] == 0


def test_summary_counts_patients_not_annotations():
    # two annotations on one patient is one patient
    found = pd.DataFrame([
        annotation("redaction", "general", entity="e1"),
        annotation("redaction", "general", entity="e2"),
    ])
    summary = summarize_annotations("Test", found, exclusion_flags(found))
    assert summary["Exclusion-grade annotations"] == 2
    assert summary["Patients with an exclusion-grade annotation"] == 1
