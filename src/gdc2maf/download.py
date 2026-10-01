"""Download the selected MAFs with gdc-client and verify them.

Only open-access files are handled, so no GDC token is needed. Every file is
checked by size and md5 against the GDC metadata both before and after the
download, so a rerun downloads only what is missing or corrupt and the check
table is the record of what the analysis actually read.
"""

import hashlib
import logging
import os
import subprocess

import pandas as pd

logger = logging.getLogger(__name__)

MANIFEST_COLUMNS = {
    "file_id": "id",
    "file_name": "filename",
    "md5sum": "md5",
    "file_size": "size",
    "state": "state",
}


def write_manifest(files, manifest_path):
    """Write a gdc-client manifest (``id filename md5 size state``) for ``files``.

    Parameters
    ----------
    files : pd.DataFrame
        Rows with ``file_id``, ``file_name``, ``md5sum``, ``file_size`` and
        ``state`` columns, e.g. the selected files from step 5.
    manifest_path : str
        Output path of the tab-separated manifest.
    """
    manifest = files[list(MANIFEST_COLUMNS)].rename(columns=MANIFEST_COLUMNS)
    manifest.to_csv(manifest_path, sep="\t", index=False)


def file_md5(path, chunk_size=1 << 20):
    """Return the hex md5 checksum of the file at ``path``."""
    md5 = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            md5.update(chunk)
    return md5.hexdigest()


def check_downloaded_files(files, download_dir):
    """Check that every file in ``files`` is on disk with the expected size and md5.

    gdc-client stores each file as ``{download_dir}/{file_id}/{file_name}``.

    Parameters
    ----------
    files : pd.DataFrame
        Rows with ``case_id``, ``submitter_id``, ``file_id``, ``file_name``,
        ``file_size`` and ``md5sum``.
    download_dir : str
        Directory given to ``gdc-client download -d``.

    Returns
    -------
    pd.DataFrame
        One row per expected file with ``path``, ``observed_size``,
        ``observed_md5`` and ``status``: ``"ok"``, ``"missing"``,
        ``"size mismatch"`` or ``"md5 mismatch"``.
    """
    rows = []
    for _, f in files.iterrows():
        path = os.path.join(download_dir, f["file_id"], f["file_name"])
        row = dict(
            case_id=f["case_id"],
            submitter_id=f["submitter_id"],
            file_id=f["file_id"],
            file_name=f["file_name"],
            path=path,
            expected_size=f["file_size"],
            expected_md5=f["md5sum"],
            observed_size=None,
            observed_md5=None,
        )
        if not os.path.isfile(path):
            row["status"] = "missing"
        else:
            row["observed_size"] = os.path.getsize(path)
            if row["observed_size"] != f["file_size"]:
                row["status"] = "size mismatch"
            else:
                row["observed_md5"] = file_md5(path)
                row["status"] = (
                    "ok" if row["observed_md5"] == f["md5sum"] else "md5 mismatch"
                )
        rows.append(row)
    return pd.DataFrame(rows)


def run_gdc_client(files, gdc_client, download_dir, work_prefix, n_clients=8):
    """Download ``files`` with ``n_clients`` gdc-client processes in parallel.

    gdc-client downloads the files of one manifest one after another (its
    ``-n`` option only splits a single file into parallel parts). For many
    small MAFs, the files are split round-robin into ``n_clients`` manifests,
    each downloaded by its own gdc-client process. GDC annotations and related
    files are not downloaded (``--no-annotations --no-related-files``).

    Parameters
    ----------
    files : pd.DataFrame
        Files to download (manifest columns, see :func:`write_manifest`).
    gdc_client : str
        Path to the gdc-client executable.
    download_dir : str
        Target directory.
    work_prefix : str
        Path prefix for the temporary manifests and per-client logs
        (``{work_prefix}_gdc_client_{i}.log``).
    n_clients : int, optional
        Number of parallel gdc-client processes.

    Returns
    -------
    list of int
        Exit code of each gdc-client process.
    """
    n_clients = max(1, min(n_clients, len(files)))
    processes = []
    for i in range(n_clients):
        manifest_path = f"{work_prefix}_manifest_to_download_{i}.txt"
        write_manifest(files.iloc[i::n_clients], manifest_path)
        log_path = f"{work_prefix}_gdc_client_{i}.log"
        processes.append(
            (
                manifest_path,
                subprocess.Popen(
                    [
                        gdc_client,
                        "download",
                        "-m",
                        manifest_path,
                        "-d",
                        download_dir,
                        "-n",
                        "1",
                        "--no-annotations",
                        "--no-related-files",
                        "--log-file",
                        log_path,
                        "--color_off",
                    ],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                ),
            )
        )
    exit_codes = []
    for manifest_path, process in processes:
        exit_codes.append(process.wait())
        os.remove(manifest_path)
    return exit_codes


def download_files(
    cohort_name,
    selected,
    out_dir,
    download_dir,
    gdc_client="tools/gdc-client/2.3/gdc-client",
    n_clients=8,
    download=True,
):
    """Download a cohort's selected MAF files and verify them.

    1. Write ``{out_dir}/{cohort_name}_manifest.txt`` for all selected files
       (the record of what the analysis uses).
    2. Check which files are already on disk and valid (size + md5).
    3. If ``download`` is True, download only missing / invalid files with
       ``n_clients`` parallel gdc-client processes (open access, no token;
       see :func:`run_gdc_client`).
    4. Check all files again and write ``{out_dir}/{cohort_name}_download_check.tsv``.

    Parameters
    ----------
    cohort_name : str
        Cohort label, e.g. ``"Lung_MALE"``.
    selected : pd.DataFrame
        Selected files from :func:`gdc2maf.selection.select_one_file_per_case`.
    out_dir : str
        Directory for the manifest and check table.
    download_dir : str
        Where gdc-client stores the files (``{download_dir}/{file_id}/{file_name}``).
    gdc_client : str, optional
        Path to the gdc-client executable.
    n_clients : int, optional
        Parallel gdc-client processes.
    download : bool, optional
        False only checks files already on disk.

    Returns
    -------
    pd.DataFrame
        The check table (see :func:`check_downloaded_files`); ``status == "ok"``
        for every file when the download is complete.
    """
    os.makedirs(download_dir, exist_ok=True)
    write_manifest(selected, f"{out_dir}/{cohort_name}_manifest.txt")

    check = check_downloaded_files(selected, download_dir)
    to_download = selected[
        selected["file_id"].isin(check.loc[check["status"] != "ok", "file_id"])
    ]
    logger.info(
        f"{cohort_name}: {len(selected)} selected files, "
        f"{len(selected) - len(to_download)} already valid on disk"
    )

    if download and len(to_download):
        work_prefix = f"{out_dir}/{cohort_name}"
        logger.info(
            f"  downloading {len(to_download)} files with "
            f"{n_clients} gdc-client processes -> {download_dir}"
        )
        exit_codes = run_gdc_client(
            to_download, gdc_client, download_dir, work_prefix, n_clients
        )
        if any(exit_codes):
            logger.warning(
                f"  gdc-client exit codes {exit_codes}; "
                f"see {work_prefix}_gdc_client_*.log"
            )
        check = check_downloaded_files(selected, download_dir)

    check.to_csv(f"{out_dir}/{cohort_name}_download_check.tsv", sep="\t", index=False)
    counts = check["status"].value_counts().to_dict()
    logger.info(f"  check: {counts}")
    return check
