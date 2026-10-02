"""Clinical data for a cohort's cases.

Auxiliary to the cohort-to-MAF path: nothing here is needed to produce the
merged MAF, but ``demographic.sex_at_birth`` is already fetched to define a
cohort, and the demographic / diagnosis / follow-up / exposure columns are what
a cohort table is usually reported with.

Multiple records per case are collapsed: the primary diagnosis (or the first
diagnosis if none is marked primary), the last follow-up by
``days_to_follow_up``, and the first exposure record. Counts of the records
seen are kept so the collapsing is visible in the output.
"""

import logging

import pandas as pd

from .api import gdc_in, gdc_query

logger = logging.getLogger(__name__)


#: ``demographic`` fields fetched for each case.
DEMOGRAPHIC_FIELDS = [
    "sex_at_birth",
    "race",
    "ethnicity",
    "country_of_residence_at_enrollment",
    "vital_status",
    "days_to_birth",
    "days_to_death",
    "age_at_index",
    "age_is_obfuscated",
]
# demographic.country_of_residence_at_enrollment is the only populated country
# field in these projects; demographic.country_of_birth is empty for every case.
COUNTRY_FIELD = "country_of_residence_at_enrollment"
DIAGNOSIS_FIELDS = [
    "age_at_diagnosis",
    "year_of_diagnosis",
    "days_to_last_follow_up",
    "primary_diagnosis",
    "morphology",
    "tissue_or_organ_of_origin",
    "ajcc_pathologic_stage",
    "ajcc_pathologic_t",
    "ajcc_pathologic_n",
    "ajcc_pathologic_m",
    "tumor_grade",
    "prior_malignancy",
    "synchronous_malignancy",
]
EXPOSURE_FIELDS = [
    "tobacco_smoking_status",
    "pack_years_smoked",
    "tobacco_smoking_onset_year",
    "tobacco_smoking_quit_year",
    "alcohol_history",
]


def gdc_clinical(case_ids):
    """Fetch one row of clinical data per case from the GDC.

    A GDC case can hold several ``diagnoses`` (primary, recurrence,
    metastasis, prior primary...) and several ``follow_ups``. This function
    keeps the diagnosis flagged ``diagnosis_is_primary_disease`` and takes the
    latest ``days_to_follow_up`` across follow-up records.

    Parameters
    ----------
    case_ids : list of str
        GDC case UUIDs.

    Returns
    -------
    pd.DataFrame
        One row per case. Day counts are relative to diagnosis, as in the GDC.
        Derived columns:

        - ``age_at_diagnosis_years``: ``floor(age_at_diagnosis / 365.25)``.
        - ``year_of_birth_est``: ``year_of_diagnosis - age_at_diagnosis_years``
          (the GDC no longer releases ``year_of_birth``; this is +/- 1 year).
        - ``days_to_last_follow_up_max``: max of ``follow_ups.days_to_follow_up``.
        - ``n_primary_diagnoses``, ``n_exposures``: record counts, for QC.
    """
    fields = (
        ["case_id", "submitter_id"]
        + [f"demographic.{f}" for f in DEMOGRAPHIC_FIELDS]
        + [f"diagnoses.{f}" for f in DIAGNOSIS_FIELDS]
        + ["diagnoses.diagnosis_is_primary_disease", "follow_ups.days_to_follow_up"]
        + [f"exposures.{f}" for f in EXPOSURE_FIELDS]
    )
    hits = gdc_query("cases", gdc_in("case_id", case_ids), fields)

    rows = []
    for h in hits:
        demographic = h.get("demographic") or {}
        primary_dx = [
            d for d in h.get("diagnoses", []) if d.get("diagnosis_is_primary_disease")
        ]
        diagnosis = primary_dx[0] if primary_dx else {}
        exposures = h.get("exposures", [])
        exposure = exposures[0] if exposures else {}
        follow_up_days = [
            fu["days_to_follow_up"]
            for fu in h.get("follow_ups", [])
            if fu.get("days_to_follow_up") is not None
        ]

        row = dict(case_id=h["case_id"], submitter_id=h["submitter_id"])
        row.update({f: demographic.get(f) for f in DEMOGRAPHIC_FIELDS})
        row.update({f: diagnosis.get(f) for f in DIAGNOSIS_FIELDS})
        row.update({f: exposure.get(f) for f in EXPOSURE_FIELDS})
        row["days_to_last_follow_up_max"] = (
            max(follow_up_days) if follow_up_days else None
        )
        row["n_primary_diagnoses"] = len(primary_dx)
        row["n_exposures"] = len(exposures)
        rows.append(row)

    clinical = pd.DataFrame(rows)
    clinical["age_at_diagnosis_years"] = (
        clinical["age_at_diagnosis"] // 365.25
    ).astype("Int64")
    clinical["year_of_birth_est"] = (
        clinical["year_of_diagnosis"] - clinical["age_at_diagnosis_years"]
    ).astype("Int64")
    return clinical


def gdc_case_country(case_ids):
    """Fetch the country of residence at enrollment for ``case_ids``.

    A small, single-field version of :func:`gdc_clinical`, so country can be
    added to clinical tables written before ``COUNTRY_FIELD`` was part of
    :data:`DEMOGRAPHIC_FIELDS` without re-running the whole clinical step.

    Parameters
    ----------
    case_ids : list of str
        GDC case UUIDs.

    Returns
    -------
    pd.DataFrame
        ``case_id``, ``submitter_id`` and ``country_of_residence_at_enrollment``
        (missing where the GDC has no value).
    """
    fields = ["case_id", "submitter_id", f"demographic.{COUNTRY_FIELD}"]
    hits = gdc_query("cases", gdc_in("case_id", case_ids), fields)
    return pd.DataFrame(
        [
            {
                "case_id": h["case_id"],
                "submitter_id": h["submitter_id"],
                COUNTRY_FIELD: (h.get("demographic") or {}).get(COUNTRY_FIELD),
            }
            for h in hits
        ]
    )


def fetch_clinical(cohort_name, cases, out_dir):
    """Fetch clinical data for ``cases`` and write it as a TSV.

    Parameters
    ----------
    cohort_name : str
        Cohort label used in the output file name, e.g. ``"Lung_MALE"``.
    cases : pd.DataFrame
        Cases from :func:`gdc2maf.cases.fetch_cases`; every ``case_id`` is
        queried, included or not, so the table can be reported alongside the
        attrition.
    out_dir : str
        Directory for ``{cohort_name}_clinical.tsv``.

    Returns
    -------
    pd.DataFrame
        One row of clinical data per case.
    """
    clinical = gdc_clinical(cases["case_id"].tolist())
    clinical = clinical.sort_values("submitter_id").reset_index(drop=True)

    out_path = f"{out_dir}/{cohort_name}_clinical.tsv"
    clinical.to_csv(out_path, sep="\t", index=False)

    logger.info(
        f"{cohort_name}: clinical for {len(clinical)} / "
        f"{len(cases)} cases -> {out_path}"
    )
    logger.info(
        f"  no primary diagnosis record: {(clinical['n_primary_diagnoses'] == 0).sum()}"
    )
    logger.info(
        f"  >1 exposure record (first kept): {(clinical['n_exposures'] > 1).sum()}"
    )
    logger.info(f"  age obfuscated: {clinical['age_is_obfuscated'].eq(True).sum()}")
    return clinical
