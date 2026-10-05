"""GDC annotations on a cohort's cases, aliquots and files.

The GDC annotates biospecimens and files with notes from the submitting centre
and from its own curation: prior treatment, a compromised sample, a genotype
mismatch, an item the centre has flagged do-not-use, and so on. A cohort should
know whether any of its selected files carries one.

**Having an annotation does not mean a file is unusable.** Every open WXS
ensemble MAF in the TCGA projects carries one ``Notification / General``
annotation reading "Variants from SomaticSniper are not included." -- a note
about how the file was made. Excluding on the presence of an annotation would
drop the entire cohort. What matters is the classification and the category,
which is why :data:`EXCLUSION_CLASSIFICATIONS` and
:data:`EXCLUSION_CATEGORIES` are narrow and
:func:`exclusion_flags` is the only thing that reads them.

Nothing here excludes anything on its own. It fetches and classifies, and the
caller decides -- the same way the merge QC flags a sample for review rather
than dropping it.
"""

import logging

import pandas as pd

from .api import gdc_in, gdc_query

logger = logging.getLogger(__name__)

#: Fields fetched for each annotation.
ANNOTATION_FIELDS = [
    "annotation_id",
    "entity_type",
    "entity_id",
    "entity_submitter_id",
    "category",
    "classification",
    "status",
    "notes",
    "created_datetime",
    "case_id",
    "case_submitter_id",
    "project.project_id",
]

#: Classifications that mean the GDC itself considers the item withdrawn.
#: ``redaction`` is the strongest signal the API carries; ``notification``,
#: ``centernotification`` and ``observation`` are commentary.
EXCLUSION_CLASSIFICATIONS = ("redaction",)

#: Categories that question the identity or integrity of the material, rather
#: than describing the patient's history or how a file was processed. Matched
#: case-insensitively. Deliberately narrow: categories such as
#: ``prior malignancy``, ``neoadjuvant therapy`` and ``item is noncanonical``
#: are clinical or curatorial facts a study may well want to keep, so they are
#: reported but not listed here.
EXCLUSION_CATEGORIES = (
    "item flagged dnu",
    "biospecimen identity unknown",
    "genotype mismatch",
    "sample compromised",
    "duplicate item",
)


def fetch_annotations(entity_ids, fields=None):
    """Fetch every GDC annotation attached to the given entities.

    Parameters
    ----------
    entity_ids : iterable of str
        UUIDs of cases, samples, aliquots, portions, analytes or files. The
        ``/annotations`` endpoint matches them against ``entity_id``, so an
        annotation on a case is only returned for the case's own UUID, not for
        its aliquots, and the other way round -- pass every level you care
        about.
    fields : list of str or None, optional
        Fields to fetch; ``None`` uses :data:`ANNOTATION_FIELDS`.

    Returns
    -------
    pd.DataFrame
        One row per annotation, with ``classification`` and ``category``
        lower-cased so they can be compared. Empty with the expected columns
        when nothing is annotated.
    """
    entity_ids = [e for e in dict.fromkeys(entity_ids) if pd.notna(e)]
    fields = list(fields or ANNOTATION_FIELDS)
    if not entity_ids:
        return pd.DataFrame(columns=fields)

    hits = gdc_query("annotations", gdc_in("entity_id", entity_ids), fields)
    rows = [
        {
            "annotation_id": h.get("annotation_id"),
            "entity_type": h.get("entity_type"),
            "entity_id": h.get("entity_id"),
            "entity_submitter_id": h.get("entity_submitter_id"),
            "classification": (h.get("classification") or "").lower(),
            "category": (h.get("category") or "").lower(),
            "status": h.get("status"),
            "notes": h.get("notes"),
            "created_datetime": h.get("created_datetime"),
            "case_id": h.get("case_id"),
            "case_submitter_id": h.get("case_submitter_id"),
            "project_id": (h.get("project") or {}).get("project_id"),
        }
        for h in hits
    ]
    return pd.DataFrame(rows)


def exclusion_flags(annotations, classifications=EXCLUSION_CLASSIFICATIONS,
                    categories=EXCLUSION_CATEGORIES):
    """The annotations that question the material, not just describe it.

    Parameters
    ----------
    annotations : pd.DataFrame
        Output of :func:`fetch_annotations`.
    classifications : iterable of str, optional
        Classifications treated as exclusion-grade; see
        :data:`EXCLUSION_CLASSIFICATIONS`.
    categories : iterable of str, optional
        Categories treated as exclusion-grade; see
        :data:`EXCLUSION_CATEGORIES`.

    Returns
    -------
    pd.DataFrame
        The subset of ``annotations`` matching either set.
    """
    if annotations.empty:
        return annotations
    classifications = {c.lower() for c in classifications}
    categories = {c.lower() for c in categories}
    matches = annotations["classification"].isin(classifications) | annotations[
        "category"
    ].isin(categories)
    return annotations[matches]


def summarize_annotations(cohort_name, annotations, flagged):
    """Summarise a cohort's annotations for a cross-cohort report.

    Parameters
    ----------
    cohort_name : str
        Cohort label.
    annotations : pd.DataFrame
        Every annotation found, from :func:`fetch_annotations`.
    flagged : pd.DataFrame
        The exclusion-grade subset, from :func:`exclusion_flags`.

    Returns
    -------
    dict
        ``Cohort``, the annotation count, the count per classification, and
        the exclusion-grade counts with the patients behind them.
    """
    summary = {
        "Cohort": cohort_name,
        "Annotations found": len(annotations),
        "Exclusion-grade annotations": len(flagged),
        "Patients with an exclusion-grade annotation": (
            flagged["case_submitter_id"].nunique() if len(flagged) else 0
        ),
    }
    if len(annotations):
        for classification, n in annotations["classification"].value_counts().items():
            summary[f"  {classification}"] = int(n)
    if len(flagged):
        for category, n in flagged["category"].value_counts().items():
            summary[f"  flagged: {category}"] = int(n)
    return summary
