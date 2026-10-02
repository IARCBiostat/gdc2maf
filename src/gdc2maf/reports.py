"""Cohort reports: attrition, and one-row-per-cohort summaries.

Each ``summarize_*`` function returns a plain dict for one cohort, so several
cohorts' summaries can be stacked side by side with :func:`summary_table`.
:func:`summarize_attrition` is the exception: it returns a table of its own,
one row per (step, reason).
"""

import logging

import pandas as pd

from .attrition import (
    STEP_DOWNLOAD,
    STEP_MAF_FILES,
    STEP_MERGE,
    STEP_PROJECT,
    STEP_SEX,
    exclude_cases,
    maf_step,
)
from .spec import WXS_ENSEMBLE_MAF
from .files import DUPLICATE_KINDS
from .selection import SELECTION_KEYS

logger = logging.getLogger(__name__)


def summarize_maf_files(cohort_name, cases, maf_files, spec=WXS_ENSEMBLE_MAF):
    """Summarise a cohort's MAF file listing for a cross-cohort report.

    Parameters
    ----------
    cohort_name : str
        Cohort label, e.g. ``"Lung_MALE"``.
    cases : pd.DataFrame
        Output of :func:`gdc2maf.cases.fetch_cases`.
    maf_files : pd.DataFrame
        Output of :func:`gdc2maf.files.list_maf_files`.

    Returns
    -------
    dict
        Summary statistics for one cohort. The analyte is the last character
        of the aliquot barcode's 5th field (``D`` native DNA; ``W``/``X``
        whole-genome amplified).
    """
    files_per_case = maf_files.groupby("case_id").size()
    tumor_analyte = maf_files["tumor_aliquot_barcode"].str.split("-").str[4].str[-1]
    # cases meeting the cohort definition (project + sex), with or without a MAF
    n_cohort = int(
        (
            cases["exclusion_step"].isna()
            | (cases["exclusion_step"] == maf_step(spec))
        ).sum()
    )

    summary = {
        "Cohort": cohort_name,
        "Cases in cohort": n_cohort,
        "Cases with a MAF": maf_files["case_id"].nunique(),
        "% cases with a MAF": round(100 * maf_files["case_id"].nunique() / n_cohort, 1),
        "MAF files": len(maf_files),
        "Patients with >1 file": int((files_per_case > 1).sum()),
        "Max files per patient": int(files_per_case.max()),
        "Distinct tumour samples": maf_files["tumor_sample_barcode"].nunique(),
        "Tumour aliquots WGA (W/X)": int(tumor_analyte.isin(["W", "X"]).sum()),
        "Total size (MB)": round(maf_files["file_size"].sum() / 1e6, 1),
    }
    for col, label in [
        ("tumor_sample_type", "Tumour"),
        ("normal_sample_type", "Normal"),
    ]:
        for sample_type, n in maf_files[col].value_counts().items():
            summary[f"{label}: {sample_type}"] = int(n)
    return summary


def summarize_selection(cohort_name, selected, dropped):
    """Summarise the one-file-per-patient selection for a cross-cohort report.

    Parameters
    ----------
    cohort_name : str
        Cohort label, e.g. ``"Lung_MALE"``.
    selected, dropped : pd.DataFrame
        Outputs of :func:`gdc2maf.selection.select_one_file_per_case`.

    Returns
    -------
    dict
        Counts of selected / dropped files, drop reasons and compromises in
        the selected files (non-primary, FFPE, WGA, solid-tissue normal).
    """
    summary = {
        "Cohort": cohort_name,
        "Patients (= selected files)": len(selected),
        "Dropped files": len(dropped),
    }
    for _, _, label in SELECTION_KEYS:
        summary[f"Dropped: {label}"] = int((dropped["drop_reason"] == label).sum())

    summary["Selected: non-primary tumour"] = int(selected["rank_primary"].sum())
    summary["Selected: FFPE tumour"] = int(selected["rank_ffpe"].sum())
    summary["Selected: WGA tumour DNA"] = int(selected["rank_native"].sum())
    summary["Selected: solid-tissue normal"] = int((selected["rank_normal"] == 1).sum())
    return summary


def summarize_download(cohort_name, selected, check):
    """Summarise the download check for a cross-cohort report.

    Parameters
    ----------
    cohort_name : str
        Cohort label, e.g. ``"Lung_MALE"``.
    selected : pd.DataFrame
        Selected files (step 5).
    check : pd.DataFrame
        Output of :func:`gdc2maf.download.download_files`.

    Returns
    -------
    dict
        Selected vs verified file counts, counts per failure status, and
        whether the verified files match the selection exactly.
    """
    ok = check[check["status"] == "ok"]
    summary = {
        "Cohort": cohort_name,
        "Selected files": len(selected),
        "Files in check": len(check),
        "Verified (size + md5)": len(ok),
        "Missing": int((check["status"] == "missing").sum()),
        "Size mismatch": int((check["status"] == "size mismatch").sum()),
        "md5 mismatch": int((check["status"] == "md5 mismatch").sum()),
        "Total size (MB)": round(ok["expected_size"].sum() / 1e6, 1),
        "Verified = selection": "yes"
        if set(ok["file_id"]) == set(selected["file_id"])
        else "NO",
    }
    return summary


def summarize_duplicates(cohort_name, duplicates):
    """Count patients per duplicate kind for a cross-cohort report.

    Parameters
    ----------
    cohort_name : str
        Cohort label, e.g. ``"Lung_MALE"``.
    duplicates : pd.DataFrame
        Output of :func:`gdc2maf.files.inspect_duplicates`.

    Returns
    -------
    dict
        ``Cohort``, number of patients with >1 file, number of files they
        hold, and the number of patients per duplicate kind.
    """
    summary = {
        "Cohort": cohort_name,
        "Patients with >1 file": len(duplicates),
        "Files for those patients": int(duplicates["n_files"].sum()),
    }
    for kind, label in DUPLICATE_KINDS.items():
        summary[label] = int(duplicates[kind].sum())
    return summary


def summarize_attrition(
    cohort_name,
    cases,
    maf_files,
    download_check=None,
    sample_qc=None,
    extra_steps=(),
    spec=WXS_ENSEMBLE_MAF,
):
    """Count the cases lost at each pipeline step, by reason.

    Parameters
    ----------
    cohort_name : str
        Cohort label, e.g. ``"Lung_MALE"``.
    cases : pd.DataFrame
        Output of :func:`gdc2maf.cases.fetch_cases` (all project cases, with
        exclusion columns).
    maf_files : pd.DataFrame
        Output of :func:`gdc2maf.files.list_maf_files`. Included cases with no file
        listed are counted as lost at the MAF file listing step.
    download_check : pd.DataFrame or None, optional
        Output of :func:`gdc2maf.download.download_files`. Cases
        whose selected file is not ``"ok"`` are counted as lost at the
        download step. ``None`` skips that step.
    sample_qc : pd.DataFrame or None, optional
        Output of :func:`gdc2maf.maf.merge_mafs`. Cases
        whose MAF has no variants are counted as lost at the merge step.
        ``None`` skips that step.
    extra_steps : sequence of (str, iterable, str), optional
        Further attrition steps to append after the merge, as
        ``(step_label, case_ids_lost, reason)`` triples. Downstream analyses
        (signature extraction, for example) live outside this package but can
        still be reported in the same table this way.

    Returns
    -------
    pd.DataFrame
        One row per (step, reason) with ``n_removed`` and ``n_remaining``,
        starting from all cases of the project(s). Steps where no case is lost
        are shown with ``n_removed = 0``.
    """
    cases = cases.copy()
    cases = exclude_cases(
        cases,
        ~cases["case_id"].isin(maf_files["case_id"]),
        STEP_MAF_FILES,
        "no MAF file returned by /files",
    )
    steps = [STEP_SEX, maf_step(spec), STEP_MAF_FILES]
    if download_check is not None:
        failed = download_check[download_check["status"] != "ok"].set_index("case_id")[
            "status"
        ]
        cases = exclude_cases(
            cases,
            cases["case_id"].isin(failed.index),
            STEP_DOWNLOAD,
            "selected file " + cases["case_id"].map(failed).fillna(""),
        )
        steps.append(STEP_DOWNLOAD)
    if sample_qc is not None:
        no_variants = sample_qc[~sample_qc["has_variants"].astype(bool)].set_index(
            "case_id"
        )
        reason = no_variants["no_variants_reason"].map(
            {
                "empty GDC MAF": "MAF has no variants (after GDC masking)",
                "only variants with removed GDC_FILTER flags": (
                    "only variants with removed GDC_FILTER flags"
                ),
            }
        )
        cases = exclude_cases(
            cases,
            cases["case_id"].isin(no_variants.index),
            STEP_MERGE,
            cases["case_id"].map(reason).fillna(""),
        )
        steps.append(STEP_MERGE)
    for step_label, case_ids_lost, reason in extra_steps:
        cases = exclude_cases(
            cases, cases["case_id"].isin(list(case_ids_lost)), step_label, reason
        )
        steps.append(step_label)

    rows = [
        dict(
            cohort=cohort_name,
            step=STEP_PROJECT,
            reason="all cases",
            n_removed=0,
            n_remaining=len(cases),
        )
    ]
    n_remaining = len(cases)
    for step in steps:
        lost = cases[cases["exclusion_step"] == step]
        if lost.empty:
            rows.append(
                dict(
                    cohort=cohort_name,
                    step=step,
                    reason="-",
                    n_removed=0,
                    n_remaining=n_remaining,
                )
            )
            continue
        for reason, n in lost["exclusion_reason"].value_counts().items():
            n_remaining -= n
            rows.append(
                dict(
                    cohort=cohort_name,
                    step=step,
                    reason=reason,
                    n_removed=int(n),
                    n_remaining=n_remaining,
                )
            )
    return pd.DataFrame(rows)


def summary_table(summaries):
    """Turn per-cohort summary dicts into one table for printing and saving.

    Parameters
    ----------
    summaries : list of dict
        Dicts with a ``"Cohort"`` key, e.g. from
        :func:`summarize_maf_files`.

    Returns
    -------
    pd.DataFrame
        One column per cohort, one row per statistic. Statistics missing for a
        cohort are 0; values are formatted as text (``1085``, ``88.3``).
    """
    return (
        pd.DataFrame(
            {
                s["Cohort"]: {k: v for k, v in s.items() if k != "Cohort"}
                for s in summaries
            }
        )
        .fillna(0)
        .map(lambda v: f"{v:g}" if isinstance(v, (int, float)) else v)
    )
