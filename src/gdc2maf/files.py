"""Listing a cohort's MAF files, and describing patients that have several.

Nothing is downloaded here. A patient can have more than one MAF file (several
tumour aliquots, several normals, several sequencing centres);
:func:`inspect_duplicates` says how those files differ, and
:mod:`gdc2maf.selection` picks one per patient.
"""

import logging

import pandas as pd

from .api import maf_file_metadata
from .spec import WXS_ENSEMBLE_MAF
from .tcga import parse_tcga_aliquot_barcode

logger = logging.getLogger(__name__)


def list_maf_files(cohort_name, cases, out_dir, spec=WXS_ENSEMBLE_MAF):
    """List the MAF files of a cohort's included cases (no download).

    Parameters
    ----------
    cohort_name : str
        Cohort label used in file names, e.g. ``"Lung_MALE"``.
    cases : pd.DataFrame
        Cases from :func:`gdc2maf.cases.fetch_cases`; only rows with
        ``included`` are queried.
    out_dir : str
        Directory for ``{cohort_name}_files_metadata.tsv``.
    spec : gdc2maf.spec.FileSpec, optional
        File kind to list; use the same spec as
        :func:`gdc2maf.cases.fetch_cases`.

    Returns
    -------
    pd.DataFrame
        One row per MAF file; see :func:`gdc2maf.api.maf_file_metadata` for the
        columns.
    """
    case_ids = cases.loc[cases["included"], "case_id"].tolist()

    files = maf_file_metadata(spec, case_ids)
    files = files.sort_values(["submitter_id", "tumor_aliquot_barcode"]).reset_index(
        drop=True
    )

    out_path = f"{out_dir}/{cohort_name}_files_metadata.tsv"
    files.to_csv(out_path, sep="\t", index=False)

    files_per_case = files.groupby("case_id").size()
    logger.info(
        f"{cohort_name}: {len(files)} MAF files for {files['case_id'].nunique()} / "
        f"{len(case_ids)} cases -> {out_path}"
    )
    logger.info(f"  cases with >1 file: {(files_per_case > 1).sum()}")
    logger.info(f"  total size: {files['file_size'].sum() / 1e6:.1f} MB")
    return files


DUPLICATE_KINDS = {
    "tumour_sample_type_differs": (
        "Tumour sample types differ (e.g. primary vs metastatic)"
    ),
    "tumour_sample_differs": "Different tumour samples (sample type or vial)",
    "same_tumour_sample_diff_aliquot": (
        "Same tumour sample, different aliquots (portion/plate)"
    ),
    "tumour_wga_and_native": "Tumour DNA: both WGA and native",
    "tumour_all_wga": "Tumour DNA: WGA only (no native option)",
    "tumour_ffpe_and_frozen": "Tumour: both FFPE and non-FFPE",
    "normal_differs": "Different normal aliquots",
    "normal_type_differs": "Normal types differ (blood vs solid tissue)",
    "same_tumour_aliquot_diff_normal": "Same tumour aliquot, different normals",
    "center_differs": "Sequencing centres differ",
}


def inspect_duplicates(cohort_name, maf_files, out_dir):
    """Describe how the MAF files of patients with more than one file differ.

    For every patient with >1 file, flags which of :data:`DUPLICATE_KINDS`
    apply. Flags are not mutually exclusive: a patient can have, e.g., two
    tumour vials and two normals. The combination of flags is stored as
    ``pattern``.

    Parameters
    ----------
    cohort_name : str
        Cohort label used in file names, e.g. ``"Lung_MALE"``.
    maf_files : pd.DataFrame
        Output of :func:`list_maf_files`.
    out_dir : str
        Directory for ``{cohort_name}_duplicate_patients.tsv``.

    Returns
    -------
    pd.DataFrame
        One row per patient with >1 file: number of files, one boolean column
        per duplicate kind, ``pattern``, and the tumour / normal aliquot
        barcodes joined with ``;``.

    Raises
    ------
    ValueError
        If two files share the same tumour and normal aliquot.
    """
    files = maf_files.copy()
    tumor = parse_tcga_aliquot_barcode(files["tumor_aliquot_barcode"])
    files["tumor_analyte"] = tumor["analyte"]
    files["center"] = tumor["center"]
    files["tumor_is_wga"] = tumor["analyte"].isin(["W", "X"])
    files["tumor_is_ffpe"] = files["tumor_sample_type"].str.contains("FFPE")

    pairs = files[["tumor_aliquot_id", "normal_aliquot_id"]]
    if pairs.duplicated().any():
        raise ValueError(f"{cohort_name}: files with identical tumour-normal pairs")

    rows = []
    for submitter_id, f in files.groupby("submitter_id"):
        if len(f) < 2:
            continue
        # does any tumour sample / aliquot appear in more than one file?
        n_per_tumour_sample = f.groupby("tumor_sample_barcode")[
            "tumor_aliquot_id"
        ].nunique()
        n_per_tumour_aliquot = f.groupby("tumor_aliquot_id")[
            "normal_aliquot_id"
        ].nunique()

        flags = {
            "tumour_sample_type_differs": f["tumor_sample_type"].nunique() > 1,
            "tumour_sample_differs": f["tumor_sample_barcode"].nunique() > 1,
            "same_tumour_sample_diff_aliquot": (n_per_tumour_sample > 1).any(),
            "tumour_wga_and_native": f["tumor_is_wga"].nunique() > 1,
            "tumour_all_wga": f["tumor_is_wga"].all(),
            "tumour_ffpe_and_frozen": f["tumor_is_ffpe"].nunique() > 1,
            "normal_differs": f["normal_aliquot_id"].nunique() > 1,
            "normal_type_differs": f["normal_sample_type"].nunique() > 1,
            "same_tumour_aliquot_diff_normal": (n_per_tumour_aliquot > 1).any(),
            "center_differs": f["center"].nunique() > 1,
        }
        rows.append(
            dict(
                submitter_id=submitter_id,
                n_files=len(f),
                **flags,
                pattern=" + ".join(k for k, v in flags.items() if v),
                tumor_aliquot_barcodes=";".join(sorted(f["tumor_aliquot_barcode"])),
                normal_aliquot_barcodes=";".join(sorted(f["normal_aliquot_barcode"])),
            )
        )

    duplicates = pd.DataFrame(
        rows,
        columns=[
            "submitter_id",
            "n_files",
            *DUPLICATE_KINDS,
            "pattern",
            "tumor_aliquot_barcodes",
            "normal_aliquot_barcodes",
        ],
    )
    out_path = f"{out_dir}/{cohort_name}_duplicate_patients.tsv"
    duplicates.to_csv(out_path, sep="\t", index=False)
    logger.info(
        f"{cohort_name}: {len(duplicates)} patients with >1 MAF file -> {out_path}"
    )
    return duplicates
