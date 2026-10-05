"""TCGA PanCanAtlas supplemental tables, downloaded and md5-checked.

Two of the PanCanAtlas publication supplements are useful alongside a GDC
cohort:

``clinical``
    ``clinical_PANCAN_patient_with_followup.tsv``, 746 columns of curated
    clinical and follow-up data per patient. The GDC's own clinical table
    (:func:`gdc2maf.clinical.fetch_clinical`) is the better default -- it is
    pinned to the GDC release the cohort was built from -- but this one carries
    variables the GDC does not. It is also what
    :func:`gdc2maf.cases.fill_sex_from_table` expects.

``quality_annotations``
    ``merged_sample_quality_annotations.tsv``, one row per
    patient x aliquot x platform, carrying the TCGA analysis working groups'
    ``Do_not_use`` flag and their pathology exclusions.

These are publication supplements, not indexed GDC data files, so the API
publishes no checksum for them: ``/files`` returns nothing for their UUIDs.
Their md5 is therefore pinned here from the files served at the time of
writing, and checked on download, so a file that changes upstream fails loudly
instead of quietly changing an analysis. :data:`PANCAN_FILES` records the UUIDs
as published on
https://gdc.cancer.gov/about-data/publications/pancanatlas
"""

import json
import logging
import os
import urllib.request

import pandas as pd

from .download import file_md5

logger = logging.getLogger(__name__)

GDC_DATA_URL = "https://api.gdc.cancer.gov/data"
PANCANATLAS_PAGE = "https://gdc.cancer.gov/about-data/publications/pancanatlas"

#: PanCanAtlas supplements this module can fetch. Each entry is keyed by a
#: short name and holds the file name, the GDC data UUID, the md5 of the file
#: served, and a one-line description. The md5 is pinned rather than fetched:
#: these are publication supplements, and the GDC ``/files`` endpoint does not
#: index them, so there is no published checksum to compare against.
PANCAN_FILES = {
    "clinical": (
        "clinical_PANCAN_patient_with_followup.tsv",
        "0fc78496-818b-4896-bd83-52db1f533c5c",
        "ffcb35edda305dd8d615497f9214eb92",
        "curated clinical and follow-up data, one row per patient",
    ),
    "quality_annotations": (
        "merged_sample_quality_annotations.tsv",
        "1a7d7be8-675d-4e60-a105-19d4121bdebf",
        "05ddd2270fb1fb24fbdc2fe9bf7384e5",
        "AWG sample quality annotations, including the Do_not_use flag",
    ),
}

#: Column of :data:`PANCAN_FILES` ``quality_annotations`` holding the flag.
DO_NOT_USE_COLUMN = "Do_not_use"
PATIENT_COLUMN = "patient_barcode"
ALIQUOT_COLUMN = "aliquot_barcode"


def fetch_pancan_file(name, dest_dir, refresh=False, check_md5=True):
    """Download one PanCanAtlas supplement and verify its md5.

    A file already in ``dest_dir`` with the expected md5 is kept, so this is
    cheap to call on every run.

    Parameters
    ----------
    name : str
        Key of :data:`PANCAN_FILES`, ``"clinical"`` or
        ``"quality_annotations"``.
    dest_dir : str
        Directory to download into; created if missing.
    refresh : bool, optional
        Download again even if a valid copy is already there.
    check_md5 : bool, optional
        Verify against the pinned md5. ``False`` accepts whatever is served,
        which only makes sense if the upstream file is known to have changed.

    Returns
    -------
    dict
        ``path``, ``name``, ``uuid``, ``url``, ``md5``, ``expected_md5``,
        ``size`` and ``source`` (``"downloaded"`` or ``"existing"``).

    Raises
    ------
    KeyError
        If ``name`` is not in :data:`PANCAN_FILES`.
    RuntimeError
        If the downloaded file's md5 is not the pinned one.
    """
    if name not in PANCAN_FILES:
        raise KeyError(
            f"unknown PanCanAtlas file {name!r}; choose from {list(PANCAN_FILES)}"
        )
    file_name, uuid, expected_md5, description = PANCAN_FILES[name]
    url = f"{GDC_DATA_URL}/{uuid}"
    os.makedirs(dest_dir, exist_ok=True)
    path = os.path.join(dest_dir, file_name)

    source = "downloaded"
    if os.path.isfile(path) and not refresh:
        if not check_md5 or file_md5(path) == expected_md5:
            logger.info(f"PanCanAtlas {name}: using {path}")
            source = "existing"
        else:
            logger.warning(f"{path} does not match the pinned md5; downloading again")

    if source == "downloaded":
        logger.info(f"PanCanAtlas {name} ({description})")
        logger.info(f"  downloading {url} -> {path}")
        urllib.request.urlretrieve(url, path)

    observed = file_md5(path)
    if check_md5 and observed != expected_md5:
        raise RuntimeError(
            f"md5 mismatch for {file_name}: expected {expected_md5}, got {observed}. "
            f"The file served by the GDC has changed since this md5 was pinned; "
            f"check {PANCANATLAS_PAGE} and update PANCAN_FILES, or pass "
            f"check_md5=False to accept it."
        )
    return {
        "path": path,
        "name": name,
        "file_name": file_name,
        "uuid": uuid,
        "url": url,
        "md5": observed,
        "expected_md5": expected_md5,
        "size": os.path.getsize(path),
        "source": source,
    }


def fetch_pancan_clinical(dest_dir, **kwargs):
    """Download ``clinical_PANCAN_patient_with_followup.tsv``.

    See :func:`fetch_pancan_file` for the parameters and the returned record.
    The path it returns is what
    :func:`gdc2maf.cases.fill_sex_from_table` and
    ``Cohort(sex_fallback_path=...)`` expect.
    """
    return fetch_pancan_file("clinical", dest_dir, **kwargs)


def fetch_pancan_quality_annotations(dest_dir, **kwargs):
    """Download ``merged_sample_quality_annotations.tsv``.

    See :func:`fetch_pancan_file` for the parameters and the returned record.
    """
    return fetch_pancan_file("quality_annotations", dest_dir, **kwargs)


def write_pancan_record(record, path):
    """Write a :func:`fetch_pancan_file` record as JSON, for provenance."""
    with open(path, "w") as f:
        json.dump(record, f, indent=2)
    return path


def load_quality_annotations(path):
    """Read the PanCanAtlas sample quality annotations table.

    Parameters
    ----------
    path : str
        Path of ``merged_sample_quality_annotations.tsv``.

    Returns
    -------
    pd.DataFrame
        One row per patient x aliquot x platform, with ``Do_not_use`` as a
        boolean. The file writes it as the strings ``True``/``False``.
    """
    annotations = pd.read_csv(path, sep="\t", low_memory=False)
    annotations[DO_NOT_USE_COLUMN] = (
        annotations[DO_NOT_USE_COLUMN].astype(str).str.strip().str.upper() == "TRUE"
    )
    return annotations


def do_not_use_patients(annotations, platform=None):
    """Patients with any ``Do_not_use`` row in the quality annotations.

    The flag is per patient x aliquot x platform, and most of the platforms are
    not exome sequencing, so a patient flagged here is not necessarily a
    patient whose WXS data is bad: the reasons include prior treatment, other
    platforms' QC, duplicates and FFPE. This returns the patients with *any*
    such row, which is the conservative reading the TCGA analysis working
    groups applied; narrow it with ``platform`` to flag on one platform only.

    Parameters
    ----------
    annotations : pd.DataFrame
        Output of :func:`load_quality_annotations`.
    platform : str or None, optional
        Keep only rows of this ``platform``. ``None`` uses every platform.

    Returns
    -------
    set of str
        TCGA patient barcodes, e.g. ``"TCGA-05-4244"``.
    """
    flagged = annotations[annotations[DO_NOT_USE_COLUMN]]
    if platform is not None:
        flagged = flagged[flagged["platform"] == platform]
    return set(flagged[PATIENT_COLUMN].dropna())
