"""Merging a cohort's MAFs, and the QC that goes with it."""

import os

import pandas as pd
import pytest
from conftest import variant, write_maf

from gdc2maf.maf import merge_mafs, read_gdc_maf, sample_qc_metrics


def test_read_gdc_maf_separates_comments_from_variants(tmp_path):
    path = tmp_path / "a.maf.gz"
    write_maf(path, [variant("T-1", "N-1", 100)])
    comments, header, maf = read_gdc_maf(str(path))
    assert comments == ["#version gdc-1.0.0", "#filedate 20240101"]
    assert header[0] == "Hugo_Symbol"
    assert len(maf) == 1


def test_read_gdc_maf_handles_a_file_with_no_variants(tmp_path):
    path = tmp_path / "empty.maf.gz"
    write_maf(path, [])
    _, _, maf = read_gdc_maf(str(path))
    assert maf.empty


def test_merge_concatenates_and_writes_the_maf(downloaded_cohort):
    maf, sample_qc, summary = merge_mafs(
        "Test",
        downloaded_cohort["selected"],
        downloaded_cohort["download_check"],
        out_dir=downloaded_cohort["out_dir"],
    )
    assert len(maf) == 6  # two patients x three variants
    assert len(sample_qc) == 2
    assert sample_qc["has_variants"].all()
    assert summary["Cohort"] == "Test"
    assert os.path.exists(os.path.join(downloaded_cohort["out_dir"], "Test.maf"))


def test_merge_respects_maf_name(downloaded_cohort, tmp_path):
    merge_mafs(
        "Test",
        downloaded_cohort["selected"],
        downloaded_cohort["download_check"],
        out_dir=downloaded_cohort["out_dir"],
        maf_name="Test_wxs_ensemble_maf.maf",
    )
    out = os.path.join(downloaded_cohort["out_dir"], "Test_wxs_ensemble_maf.maf")
    assert os.path.exists(out)


def test_merge_refuses_a_failed_download_check(downloaded_cohort):
    check = downloaded_cohort["download_check"].copy()
    check.loc[0, "status"] = "md5 mismatch"
    with pytest.raises(ValueError, match="failed the download check"):
        merge_mafs(
            "Test",
            downloaded_cohort["selected"],
            check,
            out_dir=downloaded_cohort["out_dir"],
        )


def test_merge_refuses_the_wrong_reference_build(downloaded_cohort, tmp_path):
    selected = downloaded_cohort["selected"]
    f = selected.iloc[0]
    rows = [
        variant(f["tumor_aliquot_barcode"], f["normal_aliquot_barcode"], start)
        for start in (100, 200, 300)
    ]
    rows[0]["NCBI_Build"] = "GRCh37"
    write_maf(downloaded_cohort["download_check"].loc[0, "path"], rows)
    with pytest.raises(ValueError, match="NCBI_Build"):
        merge_mafs(
            "Test",
            selected,
            downloaded_cohort["download_check"],
            out_dir=downloaded_cohort["out_dir"],
        )


def test_merge_refuses_a_barcode_that_is_not_the_selected_aliquot(downloaded_cohort):
    selected = downloaded_cohort["selected"]
    f = selected.iloc[0]
    write_maf(
        downloaded_cohort["download_check"].loc[0, "path"],
        [variant("TCGA-ZZ-9999-01A-11D-A271-08", f["normal_aliquot_barcode"], 100)],
    )
    with pytest.raises(ValueError, match="tumour barcode"):
        merge_mafs(
            "Test",
            selected,
            downloaded_cohort["download_check"],
            out_dir=downloaded_cohort["out_dir"],
        )


def test_merge_refuses_duplicate_variants_within_a_sample(downloaded_cohort):
    selected = downloaded_cohort["selected"]
    f = selected.iloc[0]
    write_maf(
        downloaded_cohort["download_check"].loc[0, "path"],
        [variant(f["tumor_aliquot_barcode"], f["normal_aliquot_barcode"], 100)] * 2,
    )
    with pytest.raises(ValueError, match="duplicate variants"):
        merge_mafs(
            "Test",
            selected,
            downloaded_cohort["download_check"],
            out_dir=downloaded_cohort["out_dir"],
        )


def test_a_patient_with_no_variants_is_reported_not_dropped(downloaded_cohort):
    write_maf(downloaded_cohort["download_check"].loc[0, "path"], [])
    _, sample_qc, _ = merge_mafs(
        "Test",
        downloaded_cohort["selected"],
        downloaded_cohort["download_check"],
        out_dir=downloaded_cohort["out_dir"],
    )
    assert len(sample_qc) == 2
    assert (~sample_qc["has_variants"]).sum() == 1
    assert "empty GDC MAF" in set(sample_qc["no_variants_reason"].dropna())


def test_gdc_filter_variants_are_kept_by_default(downloaded_cohort):
    selected = downloaded_cohort["selected"]
    f = selected.iloc[0]
    write_maf(
        downloaded_cohort["download_check"].loc[0, "path"],
        [
            variant(
                f["tumor_aliquot_barcode"],
                f["normal_aliquot_barcode"],
                100,
                gdc_filter="NonExonic",
            )
        ],
    )
    maf, _, _ = merge_mafs(
        "Test",
        selected,
        downloaded_cohort["download_check"],
        out_dir=downloaded_cohort["out_dir"],
    )
    assert (maf["GDC_FILTER"] == "NonExonic").sum() == 1


def test_listed_gdc_filter_flags_are_removed_on_request(downloaded_cohort):
    selected = downloaded_cohort["selected"]
    f = selected.iloc[0]
    write_maf(
        downloaded_cohort["download_check"].loc[0, "path"],
        [
            variant(
                f["tumor_aliquot_barcode"],
                f["normal_aliquot_barcode"],
                100,
                gdc_filter="NonExonic",
            ),
            variant(f["tumor_aliquot_barcode"], f["normal_aliquot_barcode"], 200),
        ],
    )
    maf, sample_qc, _ = merge_mafs(
        "Test",
        selected,
        downloaded_cohort["download_check"],
        out_dir=downloaded_cohort["out_dir"],
        remove_gdc_filter_flags=["NonExonic"],
    )
    assert (maf["GDC_FILTER"] == "NonExonic").sum() == 0
    assert len(maf) == 4  # one removed from the first patient, three kept elsewhere


def test_wga_tumours_are_labelled(downloaded_cohort):
    maf, sample_qc, _ = merge_mafs(
        "Test",
        downloaded_cohort["selected"],
        downloaded_cohort["download_check"],
        out_dir=downloaded_cohort["out_dir"],
    )
    assert set(maf["tumor_analyte"]) == {"D"}
    assert not maf["tumor_is_wga"].any()
    assert not sample_qc["wga"].any()


def test_sample_metrics_count_variants_and_vaf():
    maf = pd.DataFrame(
        [
            variant("T-1", "N-1", 100, t_depth=100, t_alt_count=40),
            variant("T-1", "N-1", 200, t_depth=100, t_alt_count=5),
            variant("T-1", "N-1", 300, variant_type="INS", ref="-", alt="A"),
        ]
    ).astype(str)
    metrics = sample_qc_metrics(maf)
    assert metrics.loc["T-1", "n_variants"] == 3
    assert metrics.loc["T-1", "n_snv"] == 2
    assert metrics.loc["T-1", "n_indel"] == 1
    assert metrics.loc["T-1", "frac_vaf_lt_0.1"] == pytest.approx(1 / 3)


def test_sample_metrics_fold_substitutions_onto_pyrimidines():
    maf = pd.DataFrame(
        [
            variant("T-1", "N-1", 100, ref="C", alt="T"),
            variant("T-1", "N-1", 200, ref="G", alt="A"),  # same class as C>T
        ]
    ).astype(str)
    metrics = sample_qc_metrics(maf)
    assert metrics.loc["T-1", "frac_C>T"] == pytest.approx(1.0)
