"""Resolving a cohort specification into a table of GDC cases.

A cohort is defined by GDC project(s) plus optional case-level criteria (sex at
birth). Projects are used rather than ``primary_site`` because, for example,
``primary_site="breast"`` also returns cases of other projects biopsied in the
breast.

All cases of the project(s) are kept, with the ones failing a criterion marked
as excluded (see :mod:`gdc2maf.attrition`), so the cohort's attrition can be
reported from the full starting set.
"""

import json
import logging
import os
from datetime import date

import pandas as pd

from .api import cases_with_files, gdc_cohort_cases, gdc_data_release
from .attrition import STEP_SEX, exclude_cases, maf_step
from .spec import WXS_ENSEMBLE_MAF

logger = logging.getLogger(__name__)


def fill_sex_from_table(
    cases,
    path,
    barcode_column="bcr_patient_barcode",
    sex_column="gender",
    source="PanCan",
    encoding="latin-1",
):
    """Fill missing GDC ``sex_at_birth`` from an external clinical table.

    Only cases with no ``male`` / ``female`` value in the GDC are filled. Where
    both sources have a value, disagreements are counted and logged; the GDC
    value is kept.

    The defaults match the TCGA PanCan clinical table
    (``clinical_PANCAN_patient_with_followup.tsv``), whose ``gender`` column
    holds ``MALE`` / ``FEMALE`` keyed by ``bcr_patient_barcode``.

    Parameters
    ----------
    cases : pd.DataFrame
        Output of :func:`gdc2maf.api.gdc_cohort_cases`.
    path : str
        Path to the clinical table (TSV).
    barcode_column : str, optional
        Column holding the case submitter ID / patient barcode.
    sex_column : str, optional
        Column holding the sex value; matched case-insensitively against
        ``male`` / ``female``.
    source : str, optional
        Label recorded in ``sex_at_birth_source`` for filled cases.
    encoding : str, optional
        Encoding of the table.

    Returns
    -------
    pd.DataFrame
        ``cases`` with ``sex_at_birth`` filled where possible and
        ``sex_at_birth_source`` set to ``"GDC"``, ``source`` or missing.
    """
    table = pd.read_csv(
        path, sep="\t", encoding=encoding, usecols=[barcode_column, sex_column]
    )
    table_sex = (
        table.set_index(barcode_column)[sex_column]
        .str.lower()
        .where(lambda s: s.isin(["male", "female"]))
    )
    table_sex = cases["submitter_id"].map(table_sex)

    has_gdc = cases["sex_at_birth"].isin(["male", "female"])
    use_table = ~has_gdc & table_sex.notna()
    n_conflict = int(
        (has_gdc & table_sex.notna() & (table_sex != cases["sex_at_birth"])).sum()
    )

    cases["sex_at_birth_source"] = pd.Series(dtype=object)
    cases.loc[has_gdc, "sex_at_birth_source"] = "GDC"
    cases.loc[use_table, "sex_at_birth"] = table_sex[use_table]
    cases.loc[use_table, "sex_at_birth_source"] = source

    logger.info(
        f"  sex_at_birth: {has_gdc.sum()} from GDC, "
        f"{use_table.sum()} filled from {source}, "
        f"{(~has_gdc & ~use_table).sum()} still missing; "
        f"GDC/{source} conflicts: {n_conflict}"
    )
    return cases


def fetch_cases(
    cohort_name,
    projects,
    out_dir,
    sex_at_birth=None,
    spec=WXS_ENSEMBLE_MAF,
    sex_fallback_path=None,
    sex_fallback_source="PanCan",
):
    """Query the GDC for a cohort's cases and save them with the query record.

    Parameters
    ----------
    cohort_name : str
        Cohort label used in output file names, e.g. ``"Lung_MALE"``.
    projects : list of str
        GDC project IDs, e.g. ``["TCGA-LUAD", "TCGA-LUSC"]``.
    out_dir : str
        Directory for ``{cohort_name}_cases.tsv`` and
        ``{cohort_name}_query.json``.
    sex_at_birth : list of str or None, optional
        Values of ``demographic.sex_at_birth`` to keep. ``None`` (the default)
        keeps every case whatever its sex, including cases whose sex the GDC
        does not record.
    spec : gdc2maf.spec.FileSpec, optional
        File kind a case must have at least one of to be included.
    sex_fallback_path : str or None, optional
        Clinical table used to fill missing sex (:func:`fill_sex_from_table`).
        ``None`` uses the GDC value only.
    sex_fallback_source : str, optional
        Label for the fallback source, recorded in ``sex_at_birth_source``.

    Returns
    -------
    pd.DataFrame
        One row per case of the project(s), with ``sex_at_birth_source``, a
        ``has_maf`` flag and the exclusion columns.
    """
    cases = gdc_cohort_cases(projects=projects)
    with_maf = cases_with_files(spec, cases["case_id"].tolist())
    cases["has_maf"] = cases["case_id"].isin(with_maf)
    cases = cases.sort_values("submitter_id").reset_index(drop=True)

    logger.info(f"{cohort_name}:")
    if sex_fallback_path is not None:
        cases = fill_sex_from_table(
            cases, sex_fallback_path, source=sex_fallback_source
        )
        missing_sex = f"sex_at_birth missing in GDC and {sex_fallback_source}"
    else:
        cases["sex_at_birth_source"] = (
            cases["sex_at_birth"].notna().map({True: "GDC", False: None})
        )
        missing_sex = "sex_at_birth missing in GDC"

    cases["exclusion_step"] = pd.Series(dtype=object)
    cases["exclusion_reason"] = pd.Series(dtype=object)
    if sex_at_birth is not None:
        cases = exclude_cases(
            cases,
            ~cases["sex_at_birth"].isin(sex_at_birth),
            STEP_SEX,
            ("sex_at_birth is " + cases["sex_at_birth"]).fillna(missing_sex),
        )
    cases = exclude_cases(
        cases, ~cases["has_maf"], maf_step(spec), f"no {spec.title}"
    )
    cases["included"] = cases["exclusion_step"].isna()

    os.makedirs(out_dir, exist_ok=True)
    out_path = f"{out_dir}/{cohort_name}_cases.tsv"
    cases.to_csv(out_path, sep="\t", index=False)

    # Record the exact query so the cohort can be regenerated later
    with open(f"{out_dir}/{cohort_name}_query.json", "w") as f:
        json.dump(
            dict(
                cohort=cohort_name,
                gdc_data_release=gdc_data_release(),
                query_date=date.today().isoformat(),
                n_project_cases=len(cases),
                n_cases_included=int(cases["included"].sum()),
                projects=projects,
                sex_at_birth=sex_at_birth,
                file_spec=spec.as_record(),
                sex_fallback_path=sex_fallback_path,
                n_sex_filled_from_fallback=int(
                    (cases["sex_at_birth_source"] == sex_fallback_source).sum()
                ),
            ),
            f,
            indent=2,
        )

    logger.info(
        f"{cohort_name}: {cases['included'].sum()} / {len(cases)} project cases "
        f"included -> {out_path}"
    )
    return cases
