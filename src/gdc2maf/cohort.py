"""The cohort-to-MAF path, as one object.

:class:`Cohort` holds what defines a cohort and where its outputs live, and
exposes one method per step. Each method returns DataFrames and caches its
result to ``out_dir``, so a rerun reuses what is already there instead of
re-querying the GDC or re-downloading files. Pass ``refresh`` to recompute a
step on purpose.

    >>> from gdc2maf import Cohort
    >>> cohort = Cohort(
    ...     name="Lung_MALE",
    ...     projects=["TCGA-LUAD", "TCGA-LUSC"],
    ...     sex_at_birth=["male"],
    ...     out_dir="out/Lung_MALE",
    ...     download_dir="data/GDC/Lung_MALE",
    ... )
    >>> maf, sample_qc, _ = cohort.merged_maf()     # doctest: +SKIP

Calling a late step runs the earlier ones it needs, so the line above is the
whole pipeline. :meth:`Cohort.run` does the same and returns every table.
"""

import json
import logging
import os
from collections.abc import Sequence
from dataclasses import dataclass, field

import pandas as pd

from . import reports
from .cases import fetch_cases
from .client import ensure_gdc_client
from .clinical import fetch_clinical
from .download import download_files
from .files import inspect_duplicates, list_maf_files
from .maf import merge_mafs
from .record import download_record_text, write_if_changed
from .selection import select_one_file_per_case
from .spec import WXS_ENSEMBLE_MAF, FileSpec

logger = logging.getLogger(__name__)

#: Steps whose results are cached in ``out_dir`` and can be recomputed by name
#: with ``refresh``. Duplicate inspection and file selection are not listed:
#: they are fast, local and deterministic, so they always rerun.
CACHED_STEPS = ("cases", "clinical", "files", "download", "merge")


@dataclass
class Cohort:
    """A cohort of GDC cases and the MAF built from it.

    Parameters
    ----------
    name : str
        Cohort label, used in every output file name, e.g. ``"Lung_MALE"``.
    projects : sequence of str
        GDC project IDs, e.g. ``["TCGA-LUAD", "TCGA-LUSC"]``.
    out_dir : str
        Directory for this cohort's tables, reports and merged MAF. Created if
        missing.
    download_dir : str
        Directory gdc-client downloads into, as
        ``{download_dir}/{file_id}/{file_name}``. Keep it outside ``out_dir``
        and shared between runs so downloads are reused.
    sex_at_birth : sequence of str or None, optional
        Values of ``demographic.sex_at_birth`` to keep (``"male"``,
        ``"female"``). ``None`` (the default) keeps every case whatever its
        sex, including cases whose sex the GDC does not record.
    spec : gdc2maf.spec.FileSpec, optional
        Which GDC files are the cohort's MAFs. Defaults to
        :data:`gdc2maf.spec.WXS_ENSEMBLE_MAF`.
    sex_fallback_path : str or None, optional
        Clinical table used to fill ``sex_at_birth`` where the GDC has none
        (see :func:`gdc2maf.cases.fill_sex_from_table`). ``None`` uses the GDC
        value only.
    sex_fallback_source : str, optional
        Label recorded for cases filled from that table.
    remove_gdc_filter_flags : sequence of str, optional
        ``GDC_FILTER`` flags whose variants are dropped at the merge. Empty
        (the default) keeps every flagged variant and reports the counts.
    maf_name : str or None, optional
        File name of the merged MAF inside ``out_dir``. ``None`` uses
        ``{name}_{spec.name}.maf``. Set it to keep an existing file name when
        adopting gdc2maf in a pipeline that already has merged MAFs on disk.
    refresh : sequence of str, optional
        Steps to recompute instead of reading the cached output; names from
        :data:`CACHED_STEPS`, or ``"all"``.
    """

    name: str
    projects: Sequence[str]
    out_dir: str
    download_dir: str
    sex_at_birth: Sequence[str] | None = None
    spec: FileSpec = WXS_ENSEMBLE_MAF
    sex_fallback_path: str | None = None
    sex_fallback_source: str = "PanCan"
    remove_gdc_filter_flags: Sequence[str] = ()
    maf_name: str | None = None
    refresh: Sequence[str] = field(default_factory=tuple)

    def __post_init__(self):
        """Validate ``refresh`` and create ``out_dir``."""
        unknown = set(self.refresh) - set(CACHED_STEPS) - {"all"}
        if unknown:
            raise ValueError(
                f"unknown refresh step(s) {sorted(unknown)}; "
                f"choose from {list(CACHED_STEPS)} or 'all'"
            )
        os.makedirs(self.out_dir, exist_ok=True)
        self._cache = {}

    # ---------------------------------------------------------------- paths

    def path(self, suffix):
        """Return the path of one of this cohort's outputs, ``{name}_{suffix}``."""
        return os.path.join(self.out_dir, f"{self.name}_{suffix}")

    @property
    def cases_path(self):
        """Path of the cases table."""
        return self.path("cases.tsv")

    @property
    def clinical_path(self):
        """Path of the clinical table."""
        return self.path("clinical.tsv")

    @property
    def files_path(self):
        """Path of the MAF file metadata table."""
        return self.path("files_metadata.tsv")

    @property
    def download_check_path(self):
        """Path of the download check table."""
        return self.path("download_check.tsv")

    @property
    def maf_path(self):
        """Path of the merged MAF."""
        if self.maf_name is not None:
            return os.path.join(self.out_dir, self.maf_name)
        return self.path(f"{self.spec.name}.maf")

    @property
    def sample_qc_path(self):
        """Path of the per-sample QC table."""
        return self.path("sample_qc.tsv")

    @property
    def log_path(self):
        """Path of this cohort's run log."""
        return self.path("run.log")

    @property
    def record_path(self):
        """Path of the plain-text download record, beside the downloaded MAFs."""
        return os.path.join(self.download_dir, f"{self.name}_download_record.txt")

    # ------------------------------------------------------------- caching

    def _should_run(self, step, path):
        if "all" in self.refresh or step in self.refresh:
            return True
        if not os.path.exists(path):
            return True
        logger.info(f"{self.name}: reusing {path} (refresh with refresh=['{step}'])")
        return False

    def _cached(self, step, path, build, load=None):
        if step in self._cache:
            return self._cache[step]
        if self._should_run(step, path):
            result = build()
        else:
            result = (load or (lambda: pd.read_csv(path, sep="\t")))()
        self._cache[step] = result
        return result

    # --------------------------------------------------------------- steps

    def cases(self):
        """Cases of the cohort's project(s), with the exclusion columns.

        Returns
        -------
        pd.DataFrame
            One row per case; see :func:`gdc2maf.cases.fetch_cases`.
        """
        return self._cached(
            "cases",
            self.cases_path,
            lambda: fetch_cases(
                self.name,
                projects=list(self.projects),
                out_dir=self.out_dir,
                sex_at_birth=list(self.sex_at_birth)
                if self.sex_at_birth is not None
                else None,
                spec=self.spec,
                sex_fallback_path=self.sex_fallback_path,
                sex_fallback_source=self.sex_fallback_source,
            ),
        )

    def clinical(self):
        """Clinical data for the cohort's cases (auxiliary; not needed for the MAF).

        Returns
        -------
        pd.DataFrame
            One row per case; see :func:`gdc2maf.clinical.fetch_clinical`.
        """
        return self._cached(
            "clinical",
            self.clinical_path,
            lambda: fetch_clinical(self.name, self.cases(), self.out_dir),
        )

    def maf_files(self):
        """Every MAF file of the included cases (nothing downloaded).

        Returns
        -------
        pd.DataFrame
            One row per file; see :func:`gdc2maf.files.list_maf_files`.
        """
        return self._cached(
            "files",
            self.files_path,
            lambda: list_maf_files(
                self.name, self.cases(), self.out_dir, spec=self.spec
            ),
        )

    def duplicates(self):
        """How the files of patients with more than one MAF differ.

        Always recomputed: local, deterministic and fast.

        Returns
        -------
        pd.DataFrame
            One row per patient with >1 file; see
            :func:`gdc2maf.files.inspect_duplicates`.
        """
        if "duplicates" not in self._cache:
            self._cache["duplicates"] = inspect_duplicates(
                self.name, self.maf_files(), self.out_dir
            )
        return self._cache["duplicates"]

    def selection(self):
        """One MAF file per patient, and the files dropped with the reason.

        Always recomputed: local, deterministic and fast.

        Returns
        -------
        selected : pd.DataFrame
            The chosen file per patient.
        dropped : pd.DataFrame
            Every other file, with ``drop_reason``.
        """
        if "selection" not in self._cache:
            self._cache["selection"] = select_one_file_per_case(
                self.name, self.maf_files(), self.out_dir
            )
        return self._cache["selection"]

    def gdc_client(self, version="2.3", install_dir="tools/gdc-client", path=None):
        """Locate or install gdc-client and record which binary was used.

        Writes ``{name}_gdc_client.json`` (path, version, source, md5s) next to
        the cohort's other outputs.

        Parameters
        ----------
        version : str, optional
            Pinned version, or ``"latest"``.
        install_dir : str, optional
            Parent directory holding one folder per version.
        path : str or None, optional
            Use an existing gdc-client at this path instead of installing one.

        Returns
        -------
        dict
            The record written to ``{name}_gdc_client.json``.
        """
        info = ensure_gdc_client(
            version=version, install_dir=install_dir, gdc_client_path=path
        )
        with open(self.path("gdc_client.json"), "w") as f:
            json.dump(info, f, indent=2)
        return info

    def download(
        self,
        gdc_client=None,
        n_clients=8,
        download=True,
        version="2.3",
        install_dir="tools/gdc-client",
    ):
        """Download the selected files and verify every one by size and md5.

        Parameters
        ----------
        gdc_client : str or None, optional
            Path to a gdc-client executable. ``None`` locates or installs one
            via :meth:`gdc_client` (skipped entirely when ``download`` is
            False).
        n_clients : int, optional
            Parallel gdc-client processes.
        download : bool, optional
            False only verifies files already on disk, downloading nothing.
        version, install_dir : str, optional
            Passed to :meth:`gdc_client` when one has to be installed.

        Returns
        -------
        pd.DataFrame
            The check table; ``status == "ok"`` for every file when the
            download is complete.
        """

        def build():
            client = gdc_client
            if client is None and download:
                client = self.gdc_client(version=version, install_dir=install_dir)[
                    "path"
                ]
            selected, _ = self.selection()
            return download_files(
                self.name,
                selected,
                out_dir=self.out_dir,
                download_dir=self.download_dir,
                gdc_client=client,
                n_clients=n_clients,
                download=download,
            )

        return self._cached("download", self.download_check_path, build)

    def merged_maf(self):
        """Merge the cohort's MAFs into one MAF and run QC.

        Returns
        -------
        maf : pd.DataFrame or None
            The merged MAF, or ``None`` when a cached run is reused (the file
            is on disk at :attr:`maf_path`; read it with
            :meth:`read_merged_maf`).
        sample_qc : pd.DataFrame
            One row per selected patient, with the QC metrics and flags.
        qc_summary : dict or None
            One-row summary for :func:`gdc2maf.reports.summary_table`;
            ``None`` when a cached run is reused.
        """

        def build():
            selected, _ = self.selection()
            return merge_mafs(
                self.name,
                selected,
                self.download(),
                out_dir=self.out_dir,
                remove_gdc_filter_flags=list(self.remove_gdc_filter_flags),
                maf_name=os.path.basename(self.maf_path),
            )

        def load():
            logger.info(
                f"{self.name}: reusing {self.maf_path} and {self.sample_qc_path}"
            )
            return None, pd.read_csv(self.sample_qc_path, sep="\t"), None

        return self._cached("merge", self.sample_qc_path, build, load)

    def read_merged_maf(self, usecols=None):
        """Read the merged MAF from disk.

        Parameters
        ----------
        usecols : list of str or None, optional
            Columns to read; ``None`` reads all of them. The merged MAF of a
            large cohort is hundreds of MB, so pass the columns you need.

        Returns
        -------
        pd.DataFrame
            The merged MAF.
        """
        if usecols is None:
            return pd.read_csv(self.maf_path, sep="\t", low_memory=False)
        return pd.read_csv(
            self.maf_path,
            sep="\t",
            dtype=str,
            usecols=lambda c: c in set(usecols),
        )[list(usecols)]

    # ------------------------------------------------------------- reports

    def attrition(self, extra_steps=()):
        """Cases lost at each step, by reason, from all cases of the project(s).

        Parameters
        ----------
        extra_steps : sequence of (str, iterable, str), optional
            Further steps to append, as ``(step_label, case_ids_lost,
            reason)``; see :func:`gdc2maf.reports.summarize_attrition`.

        Returns
        -------
        pd.DataFrame
            One row per (step, reason), saved as ``{name}_case_attrition.tsv``.
        """
        _, sample_qc, _ = (
            self.merged_maf() if self._merge_available() else (None, None, None)
        )
        table = reports.summarize_attrition(
            self.name,
            self.cases(),
            self.maf_files(),
            download_check=self._download_check(),
            sample_qc=sample_qc,
            extra_steps=extra_steps,
            spec=self.spec,
        )
        table.to_csv(self.path("case_attrition.tsv"), sep="\t", index=False)
        return table

    def _merge_available(self):
        return "merge" in self._cache or os.path.exists(self.sample_qc_path)

    def _download_check(self):
        """Return the download check table, from this run or an earlier one."""
        if "download" in self._cache:
            return self._cache["download"]
        if os.path.exists(self.download_check_path):
            return pd.read_csv(self.download_check_path, sep="\t")
        return None

    def summaries(self):
        """One-row summaries of every step that has run, as a dict of dicts.

        Returns
        -------
        dict
            ``{report_name: summary_dict}``, each saved as
            ``{name}_{report_name}.tsv``. Steps that have not run are absent.
        """
        selected, dropped = self.selection()
        out = {
            "maf_files_summary": reports.summarize_maf_files(
                self.name, self.cases(), self.maf_files(), spec=self.spec
            ),
            "duplicates_summary": reports.summarize_duplicates(
                self.name, self.duplicates()
            ),
            "selection_summary": reports.summarize_selection(
                self.name, selected, dropped
            ),
        }
        download_check = self._download_check()
        if download_check is not None:
            out["download_summary"] = reports.summarize_download(
                self.name, selected, download_check
            )
        if "merge" in self._cache and self._cache["merge"][2] is not None:
            out["merge_qc_summary"] = self._cache["merge"][2]

        for report_name, summary in out.items():
            table = reports.summary_table([summary])
            table.to_csv(self.path(f"{report_name}.tsv"), sep="\t")
        return out

    def samples_for_review(self):
        """Return the samples flagged at the merge QC for manual review.

        Returns
        -------
        pd.DataFrame
            The flagged rows of the per-sample QC table, with the evidence
            columns, saved as ``{name}_samples_for_review.tsv``.
        """
        _, sample_qc, _ = self.merged_maf()
        cols = [
            "Tumor_Sample_Barcode",
            "n_variants",
            "high_count_threshold",
            "median_vaf",
            "frac_vaf_lt_0.1",
            "frac_weak_caller_support",
            "frac_C>A",
            "frac_C>G",
            "frac_C>T",
            "frac_T>C",
            "wga",
            "high_count",
            "low_vaf_excess",
            "weak_caller_support",
        ]
        cols = [c for c in cols if c in sample_qc.columns]
        flagged = sample_qc.loc[sample_qc["review"].eq(True), cols]
        flagged.to_csv(self.path("samples_for_review.tsv"), sep="\t", index=False)
        return flagged

    def write_download_record(self, record_dir=None):
        """Write the plain-text cases-and-files record next to the MAF files.

        The record names the cohort, the GDC release, the gdc-client used, how
        many cases were excluded at each step and with what reason, how many
        files were downloaded and verified, and (once the merge has run) how
        many patients ended up with variants.

        It is rewritten only when its content would change: the text holds no
        timestamp, so a rerun that downloads no new case leaves the file
        exactly as it was. When the content does change, the previous record is
        kept as a dated ``_superseded_`` copy.

        Parameters
        ----------
        record_dir : str or None, optional
            Directory for the record. ``None`` puts it in
            :attr:`download_dir`, beside the downloaded MAFs; pass
            :attr:`out_dir` to keep it with the merged MAF instead.

        Returns
        -------
        dict
            ``status`` (``"unchanged"``, ``"created"`` or ``"updated"``),
            ``path`` and ``superseded``; see
            :func:`gdc2maf.record.write_if_changed`.
        """
        selected, _ = self.selection()
        download_check = self._download_check()
        if download_check is None:
            raise ValueError(
                f"{self.name}: no download check yet; run download() first"
            )
        sample_qc = self.merged_maf()[1] if self._merge_available() else None
        text = download_record_text(
            self.name,
            self.spec,
            self.attrition(),
            selected,
            download_check,
            sample_qc=sample_qc,
            out_dir=self.out_dir,
            download_dir=self.download_dir,
        )
        path = (
            self.record_path
            if record_dir is None
            else os.path.join(record_dir, f"{self.name}_download_record.txt")
        )
        return write_if_changed(path, text)

    # ------------------------------------------------------------- the run

    def run(
        self,
        record_dir=None,
        clinical=True,
        download=True,
        merge=True,
        gdc_client=None,
        n_clients=8,
        version="2.3",
        install_dir="tools/gdc-client",
    ):
        """Run the whole path from cohort specification to merged MAF.

        Parameters
        ----------
        record_dir : str or None, optional
            Directory for the plain-text download record; ``None`` puts it
            beside the downloaded MAFs. See :meth:`write_download_record`.
        clinical : bool, optional
            Also fetch the clinical table (auxiliary to the MAF).
        download : bool, optional
            Download missing or invalid files. False only verifies what is on
            disk, which is what to use when the files are already there.
        merge : bool, optional
            Merge the verified MAFs and run QC.
        gdc_client, n_clients, version, install_dir
            Passed to :meth:`download`.

        Returns
        -------
        dict
            ``cases``, ``clinical``, ``maf_files``, ``duplicates``,
            ``selected``, ``dropped``, ``download_check``, ``sample_qc``,
            ``attrition``, ``summaries``, ``samples_for_review`` and
            ``record``, with the keys of steps that did not run set to
            ``None``.
        """
        sex = ", ".join(self.sex_at_birth) if self.sex_at_birth else "any"
        logger.info(
            f"Cohort: {self.name} (projects: {', '.join(self.projects)}; "
            f"sex at birth: {sex})"
        )
        logger.info(f"Files: {self.spec.title}")
        logger.info(f"Outputs: {self.out_dir}; downloads: {self.download_dir}")

        result = dict.fromkeys(
            [
                "cases",
                "clinical",
                "maf_files",
                "duplicates",
                "selected",
                "dropped",
                "download_check",
                "sample_qc",
                "attrition",
                "summaries",
                "samples_for_review",
                "record",
            ]
        )
        result["cases"] = self.cases()
        if clinical:
            result["clinical"] = self.clinical()
        result["maf_files"] = self.maf_files()
        result["duplicates"] = self.duplicates()
        result["selected"], result["dropped"] = self.selection()
        result["download_check"] = self.download(
            gdc_client=gdc_client,
            n_clients=n_clients,
            download=download,
            version=version,
            install_dir=install_dir,
        )
        if merge:
            _, result["sample_qc"], _ = self.merged_maf()
            result["samples_for_review"] = self.samples_for_review()
        result["summaries"] = self.summaries()
        result["attrition"] = self.attrition()
        result["record"] = self.write_download_record(record_dir=record_dir)
        return result
