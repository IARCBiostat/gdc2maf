"""Shared fixtures: synthetic GDC metadata and MAF files, no network access."""

import gzip

import pandas as pd
import pytest

from gdc2maf.selection import select_one_file_per_case

MAF_COLUMNS = [
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode",
    "Matched_Norm_Sample_Barcode",
    "NCBI_Build",
    "t_depth",
    "t_alt_count",
    "callers",
    "GDC_FILTER",
]


def maf_file_row(submitter_id, tumor_barcode, normal_barcode, file_id,
                 tumor_sample_type="Primary Tumor",
                 normal_sample_type="Blood Derived Normal"):
    """Build one row of file metadata as the GDC /files query returns it."""
    return {
        "file_id": file_id,
        "file_name": f"{file_id}.maf.gz",
        "file_size": 100,
        "md5sum": "x" * 32,
        "state": "released",
        "created_datetime": "2024-01-01T00:00:00",
        "case_id": f"case-{submitter_id}",
        "submitter_id": submitter_id,
        "project_id": "TCGA-TEST",
        "tumor_aliquot_id": f"tumor-{file_id}",
        "tumor_aliquot_barcode": tumor_barcode,
        "tumor_sample_barcode": tumor_barcode[:16],
        "tumor_sample_type": tumor_sample_type,
        "normal_aliquot_id": f"normal-{file_id}",
        "normal_aliquot_barcode": normal_barcode,
        "normal_sample_barcode": normal_barcode[:16],
        "normal_sample_type": normal_sample_type,
    }


@pytest.fixture
def maf_files():
    """Three patients: one single file, one with two aliquots, one WGA-only."""
    return pd.DataFrame(
        [
            maf_file_row(
                "TCGA-AA-0001",
                "TCGA-AA-0001-01A-11D-A271-08",
                "TCGA-AA-0001-10A-01D-A271-08",
                "f1",
            ),
            # same patient, native DNA (A271) and WGA (W) aliquots
            maf_file_row(
                "TCGA-AA-0002",
                "TCGA-AA-0002-01A-11D-A271-08",
                "TCGA-AA-0002-10A-01D-A271-08",
                "f2",
            ),
            maf_file_row(
                "TCGA-AA-0002",
                "TCGA-AA-0002-01A-11W-A271-08",
                "TCGA-AA-0002-10A-01D-A271-08",
                "f3",
            ),
            # metastatic (06) only
            maf_file_row(
                "TCGA-AA-0003",
                "TCGA-AA-0003-06A-11D-A271-08",
                "TCGA-AA-0003-10A-01D-A271-08",
                "f4",
                tumor_sample_type="Metastatic",
            ),
        ]
    )


def write_maf(path, rows):
    """Write a gzipped GDC-style MAF with comment lines and ``rows``."""
    lines = ["#version gdc-1.0.0", "#filedate 20240101", "\t".join(MAF_COLUMNS)]
    lines += ["\t".join(str(row.get(c, "")) for c in MAF_COLUMNS) for row in rows]
    with gzip.open(path, "wt") as f:
        f.write("\n".join(lines) + "\n")


def variant(tumor_barcode, normal_barcode, start, variant_type="SNP",
            ref="C", alt="T", t_depth=100, t_alt_count=40,
            callers="mutect2;muse;varscan2", gdc_filter=""):
    """Build one MAF variant row."""
    return {
        "Hugo_Symbol": "GENE",
        "Chromosome": "chr1",
        "Start_Position": start,
        "End_Position": start,
        "Variant_Type": variant_type,
        "Reference_Allele": ref,
        "Tumor_Seq_Allele2": alt,
        "Tumor_Sample_Barcode": tumor_barcode,
        "Matched_Norm_Sample_Barcode": normal_barcode,
        "NCBI_Build": "GRCh38",
        "t_depth": t_depth,
        "t_alt_count": t_alt_count,
        "callers": callers,
        "GDC_FILTER": gdc_filter,
    }


@pytest.fixture
def downloaded_cohort(tmp_path, maf_files):
    """A selection of two patients with MAFs on disk and a passing check table.

    Returns
    -------
    dict
        ``selected``, ``download_check`` and ``out_dir``, ready for
        :func:`gdc2maf.maf.merge_mafs`.
    """
    # merge_mafs works on the ranked selection table, so go through the
    # selection step rather than hand-building it
    two_patients = maf_files[maf_files["file_id"].isin(["f1", "f4"])]
    selected, _ = select_one_file_per_case(
        "Fixture", two_patients.reset_index(drop=True), str(tmp_path)
    )
    rows = []
    for _, f in selected.iterrows():
        path = tmp_path / f"{f['file_id']}.maf.gz"
        write_maf(
            path,
            [
                variant(f["tumor_aliquot_barcode"], f["normal_aliquot_barcode"], start)
                for start in (100, 200, 300)
            ],
        )
        rows.append({"file_id": f["file_id"], "path": str(path), "status": "ok"})
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    return {
        "selected": selected,
        "download_check": pd.DataFrame(rows),
        "out_dir": str(out_dir),
    }
