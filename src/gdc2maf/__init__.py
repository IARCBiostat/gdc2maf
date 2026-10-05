"""gdc2maf: from a GDC cohort specification to one verified, merged MAF.

Say which GDC project(s) and which cases you want; gdc2maf resolves the cohort,
lists the matching open-access MAF files, picks one file per patient, installs
a verified gdc-client, downloads and md5-checks every file, and merges them into
a single MAF with per-sample QC. Every patient lost on the way is accounted for
in an attrition table, with the reason.

The short version::

    from gdc2maf import Cohort, configure_logging

    configure_logging()
    cohort = Cohort(
        name="Lung_MALE",
        projects=["TCGA-LUAD", "TCGA-LUSC"],
        sex_at_birth=["male"],
        out_dir="out/Lung_MALE",
        download_dir="downloads",
    )
    result = cohort.run()
    print(result["attrition"])

or, equivalently, from a shell::

    gdc2maf maf --project TCGA-LUAD TCGA-LUSC --sex-at-birth male
                --name Lung_MALE --out out

Signature extraction, plotting and other downstream analysis are deliberately
out of scope: this package stops at the MAF.

Copyright (C) 2026 Ali Farnudi. gdc2maf is free software under the GNU
General Public License, version 3 or later, and comes with NO WARRANTY. See
the LICENSE file, or <https://www.gnu.org/licenses/gpl-3.0.html>.
"""

import logging

from .cases import fetch_cases, fill_sex_from_table
from .client import ensure_gdc_client, install_gdc_client
from .clinical import fetch_clinical
from .annotations import exclusion_flags, fetch_annotations
from .cohort import CACHED_STEPS, Cohort
from .download import check_downloaded_files, download_files, file_md5, write_manifest
from .files import inspect_duplicates, list_maf_files
from .maf import merge_mafs, read_gdc_maf
from .record import download_record_text, write_if_changed
from .pancan import (
    fetch_pancan_clinical,
    fetch_pancan_file,
    fetch_pancan_quality_annotations,
    do_not_use_patients,
    load_quality_annotations,
)
from .provenance import configure_logging, environment_info, log_environment, tee_stdout
from .selection import select_one_file_per_case
from .spec import WXS_ENSEMBLE_MAF, FileSpec

__version__ = "0.1.0"

__all__ = [
    "__version__",
    # the usual way in
    "Cohort",
    "CACHED_STEPS",
    # what counts as the cohort's MAFs
    "FileSpec",
    "WXS_ENSEMBLE_MAF",
    # the steps, for callers that want them one at a time
    "fetch_cases",
    "fill_sex_from_table",
    # GDC curation and the PanCanAtlas quality annotations
    "fetch_annotations",
    "exclusion_flags",
    "fetch_pancan_file",
    "fetch_pancan_clinical",
    "fetch_pancan_quality_annotations",
    "load_quality_annotations",
    "do_not_use_patients",
    "fetch_clinical",
    "list_maf_files",
    "inspect_duplicates",
    "select_one_file_per_case",
    "ensure_gdc_client",
    "install_gdc_client",
    "write_manifest",
    "check_downloaded_files",
    "download_files",
    "file_md5",
    "merge_mafs",
    "read_gdc_maf",
    "download_record_text",
    "write_if_changed",
    # logging and run records
    "configure_logging",
    "log_environment",
    "environment_info",
    "tee_stdout",
]

# A library should not configure logging for its caller; configure_logging()
# attaches the handlers when the caller wants to see the messages.
logging.getLogger(__name__).addHandler(logging.NullHandler())
