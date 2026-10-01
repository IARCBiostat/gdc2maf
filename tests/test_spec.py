"""The file spec that decides which GDC files are a cohort's MAFs."""

from gdc2maf import WXS_ENSEMBLE_MAF, FileSpec


def test_title_defaults_to_the_name_in_words():
    spec = FileSpec(name="wgs_ensemble_maf")
    assert spec.title == "wgs ensemble maf"


def test_explicit_title_is_kept():
    assert WXS_ENSEMBLE_MAF.title == "open WXS ensemble MAF"


def test_only_constrained_fields_become_clauses():
    spec = FileSpec(name="partial", data_type=("Masked Somatic Mutation",))
    clauses = spec.clauses()
    assert len(clauses) == 1
    assert clauses[0]["content"]["field"] == "data_type"


def test_wxs_spec_constrains_the_five_gdc_fields():
    fields = [c["content"]["field"] for c in WXS_ENSEMBLE_MAF.clauses()]
    assert fields == [
        "data_category",
        "data_type",
        "experimental_strategy",
        "analysis.workflow_type",
        "access",
    ]


def test_filter_puts_the_cases_first_and_ands_everything():
    query = WXS_ENSEMBLE_MAF.filter(["case-1", "case-2"])
    assert query["op"] == "and"
    assert query["content"][0]["content"] == {
        "field": "cases.case_id",
        "value": ["case-1", "case-2"],
    }
    assert len(query["content"]) == 6


def test_extra_conditions_are_included():
    spec = FileSpec(name="x", extra=(("platform", ("Illumina",)),))
    assert spec.clauses()[0]["content"]["field"] == "platform"


def test_as_record_is_json_serialisable():
    import json

    record = WXS_ENSEMBLE_MAF.as_record()
    assert json.loads(json.dumps(record))["name"] == "wxs_ensemble_maf"
