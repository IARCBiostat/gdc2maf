"""Tracking which cases are lost, where, and why.

Every case of the requested project(s) stays in the cases table for the whole
run. A case that fails a criterion gets ``exclusion_step`` and
``exclusion_reason`` filled in and ``included`` set to ``False``, rather than
being dropped, so that :func:`gdc2maf.reports.summarize_attrition` can report
the cohort's attrition from the full starting set.

Only a case's *first* exclusion is recorded, so each case is counted once, at
the step where it was lost.
"""

import pandas as pd

#: Attrition step labels, in pipeline order.
STEP_PROJECT = "1. GDC project(s)"
STEP_SEX = "2. sex at birth"
# Step 3 has no constant: its label names the file kind the case had to have,
# so that an attrition table says which MAF was looked for. See maf_step().
STEP_MAF_FILES = "4. MAF file listing"
STEP_DOWNLOAD = "5. download (size + md5 check)"
STEP_MERGE = "6. merged MAF"


def maf_step(spec):
    """Return step 3's label for ``spec``, e.g. ``"3. open WXS ensemble MAF"``.

    Parameters
    ----------
    spec : gdc2maf.spec.FileSpec
        The file kind a case is required to have.

    Returns
    -------
    str
        The attrition step label.
    """
    return f"3. {spec.title}"


def exclude_cases(cases, mask, step, reason):
    """Mark still-included ``cases`` matching ``mask`` as excluded.

    Only the first exclusion of a case is recorded, so each case is counted
    once, at the step where it was lost.

    Parameters
    ----------
    cases : pd.DataFrame
        Cases table with ``exclusion_step`` and ``exclusion_reason`` columns.
    mask : pd.Series of bool
        Cases to exclude.
    step : str
        Pipeline step at which the cases are lost.
    reason : str or pd.Series of str
        Why the cases are lost; a Series gives a per-case reason.

    Returns
    -------
    pd.DataFrame
        ``cases`` with the exclusion columns filled in.
    """
    mask = mask & cases["exclusion_step"].isna()
    cases.loc[mask, "exclusion_step"] = step
    cases.loc[mask, "exclusion_reason"] = (
        reason[mask] if isinstance(reason, pd.Series) else reason
    )
    return cases
