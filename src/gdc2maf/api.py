"""Thin client for the GDC REST API.

Filter helpers, a paging query function, and the cohort-level queries the
pipeline needs: a project's cases, which cases have a matching MAF, and the
metadata of those MAF files. Which files count as "the cohort's MAFs" is not
decided here; it comes in as a :class:`gdc2maf.spec.FileSpec`.
"""

import logging

import pandas as pd
import requests

logger = logging.getLogger(__name__)


GDC_API = "https://api.gdc.cancer.gov"


def gdc_in(field, values):
    """Build a GDC ``in`` filter clause for ``field`` matching any of ``values``."""
    return {"op": "in", "content": {"field": field, "value": list(values)}}


def gdc_and(*clauses):
    """Combine GDC filter clauses with a logical AND."""
    return {"op": "and", "content": list(clauses)}


def gdc_data_release():
    """Return the current GDC data release string (e.g. ``"Data Release 43.0"``)."""
    r = requests.get(f"{GDC_API}/status", timeout=60)
    r.raise_for_status()
    return r.json()["data_release"]


def gdc_query(endpoint, filters, fields, page_size=1000):
    """Query a GDC endpoint and return all hits, following pagination.

    Parameters
    ----------
    endpoint : str
        GDC endpoint name, e.g. ``"cases"`` or ``"files"``.
    filters : dict
        GDC filter object.
    fields : list of str
        Fields to return for each hit.
    page_size : int, optional
        Number of hits requested per page. Default is 1000.

    Returns
    -------
    list of dict
        All hits across pages.
    """
    hits = []
    start = 0
    while True:
        r = requests.post(
            f"{GDC_API}/{endpoint}",
            json={
                "filters": filters,
                "fields": ",".join(fields),
                "format": "JSON",
                "size": page_size,
                "from": start,
            },
            timeout=120,
        )
        r.raise_for_status()
        data = r.json()["data"]
        hits.extend(data["hits"])
        start += page_size
        if start >= data["pagination"]["total"]:
            return hits


def gdc_cohort_cases(projects, sex_at_birth=None):
    """Fetch TCGA cases for a cohort defined by project(s) and sex at birth.

    Parameters
    ----------
    projects : list of str
        GDC project IDs, e.g. ``["TCGA-LUAD", "TCGA-LUSC"]``.
    sex_at_birth : list of str or None, optional
        Values of ``demographic.sex_at_birth`` to keep (``"female"``,
        ``"male"``). Pass ``None`` to keep all cases.

    Returns
    -------
    pd.DataFrame
        One row per case with ``case_id``, ``submitter_id``, ``project_id``,
        ``primary_site``, ``disease_type`` and ``sex_at_birth`` columns.
    """
    clauses = [gdc_in("project.project_id", projects)]
    if sex_at_birth is not None:
        clauses.append(gdc_in("demographic.sex_at_birth", sex_at_birth))

    fields = [
        "case_id",
        "submitter_id",
        "project.project_id",
        "primary_site",
        "disease_type",
        "demographic.sex_at_birth",
    ]
    hits = gdc_query("cases", gdc_and(*clauses), fields)

    return pd.DataFrame(
        [
            dict(
                case_id=h["case_id"],
                submitter_id=h["submitter_id"],
                project_id=h["project"]["project_id"],
                primary_site=h.get("primary_site"),
                disease_type=h.get("disease_type"),
                sex_at_birth=(h.get("demographic") or {}).get("sex_at_birth"),
            )
            for h in hits
        ]
    )


def cases_with_files(spec, case_ids):
    """Return the subset of ``case_ids`` that have a file matching ``spec``.

    Queries ``/files`` (not ``/cases``) so that all of the spec's conditions
    are required to hold for the same file.

    Parameters
    ----------
    spec : gdc2maf.spec.FileSpec
        The file kind to look for, e.g.
        :data:`gdc2maf.spec.WXS_ENSEMBLE_MAF`.
    case_ids : list of str
        GDC case UUIDs.

    Returns
    -------
    set of str
        Case UUIDs with at least one matching file.
    """
    hits = gdc_query("files", spec.filter(case_ids), ["cases.case_id"])
    return {c["case_id"] for h in hits for c in h["cases"]}


def maf_file_metadata(spec, case_ids):
    """Fetch metadata for every file of ``case_ids`` matching ``spec``.

    Each ensemble MAF is one tumour-normal pair: it links to exactly two
    aliquots (``associated_entities``), and the sample each aliquot comes from
    carries ``tissue_type`` ``"Tumor"`` or ``"Normal"``. That GDC annotation is
    used to label the pair, rather than parsing sample-type codes from the
    barcode.

    Parameters
    ----------
    spec : gdc2maf.spec.FileSpec
        The file kind to list, e.g. :data:`gdc2maf.spec.WXS_ENSEMBLE_MAF`.
    case_ids : list of str
        GDC case UUIDs.

    Returns
    -------
    pd.DataFrame
        One row per file: file ID, name, size, md5, state, case, and the
        tumour and normal aliquot UUID, barcode, sample barcode and sample
        type.

    Raises
    ------
    ValueError
        If a file does not link to exactly one case, one tumour aliquot and
        one normal aliquot.
    """
    fields = [
        "file_id",
        "file_name",
        "file_size",
        "md5sum",
        "state",
        "created_datetime",
        "cases.case_id",
        "cases.submitter_id",
        "cases.project.project_id",
        "associated_entities.entity_id",
        "associated_entities.entity_submitter_id",
        "associated_entities.entity_type",
        "cases.samples.submitter_id",
        "cases.samples.sample_type",
        "cases.samples.tissue_type",
        "cases.samples.portions.analytes.aliquots.aliquot_id",
    ]
    hits = gdc_query("files", spec.filter(case_ids), fields)

    rows = []
    for h in hits:
        if len(h["cases"]) != 1:
            raise ValueError(f"{h['file_id']}: links to {len(h['cases'])} cases")
        case = h["cases"][0]

        # aliquot UUID -> the sample it was taken from
        aliquot_sample = {
            aliquot["aliquot_id"]: sample
            for sample in case.get("samples", [])
            for portion in sample.get("portions", [])
            for analyte in portion.get("analytes", [])
            for aliquot in analyte.get("aliquots", [])
        }

        row = dict(
            file_id=h["file_id"],
            file_name=h["file_name"],
            file_size=h["file_size"],
            md5sum=h["md5sum"],
            state=h["state"],
            created_datetime=h["created_datetime"],
            case_id=case["case_id"],
            submitter_id=case["submitter_id"],
            project_id=case["project"]["project_id"],
        )
        for entity in h.get("associated_entities", []):
            if entity["entity_type"] != "aliquot":
                continue
            sample = aliquot_sample[entity["entity_id"]]
            role = sample["tissue_type"].lower()  # "tumor" or "normal"
            if f"{role}_aliquot_id" in row:
                raise ValueError(f"{h['file_id']}: more than one {role} aliquot")
            row[f"{role}_aliquot_id"] = entity["entity_id"]
            row[f"{role}_aliquot_barcode"] = entity["entity_submitter_id"]
            row[f"{role}_sample_barcode"] = sample["submitter_id"]
            row[f"{role}_sample_type"] = sample["sample_type"]

        if "tumor_aliquot_id" not in row or "normal_aliquot_id" not in row:
            raise ValueError(f"{h['file_id']}: missing tumour or normal aliquot")
        rows.append(row)

    return pd.DataFrame(rows)
