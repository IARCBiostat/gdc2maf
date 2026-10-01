"""A plain-text record of what was downloaded, written next to the MAF files.

The downloaded MAFs are just numbered directories of ``.maf.gz`` files; nothing
in them says which cohort they belong to or how many patients were excluded on
the way. :func:`download_record_text` renders that as readable text and
:func:`write_if_changed` puts it beside the files.

The record is written **only when its content would change**. It carries no
timestamp of its own and nothing else that varies run to run, so a rerun that
downloads no new case produces exactly the same text, sees that the file already
says it, and leaves the file untouched — same bytes, same modification time.
When the content does change (a new GDC release, a wider cohort, a file that
failed verification), the previous record is kept as a dated ``_superseded_``
copy rather than lost.
"""

import json
import logging
import os
from datetime import date

logger = logging.getLogger(__name__)


def _format_size(n_bytes):
    """Return a byte count as a short human-readable string."""
    size = float(n_bytes)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size < 1024 or unit == "TB":
            return f"{size:.1f} {unit}" if unit != "B" else f"{size:.0f} B"
        size /= 1024
    return f"{size:.1f} TB"


def _query_record(out_dir, cohort_name):
    """Read ``{cohort_name}_query.json`` if the cases step wrote one."""
    path = os.path.join(out_dir, f"{cohort_name}_query.json")
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return json.load(f)


def _gdc_client_record(out_dir, cohort_name):
    """Read ``{cohort_name}_gdc_client.json`` if a download wrote one."""
    path = os.path.join(out_dir, f"{cohort_name}_gdc_client.json")
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return json.load(f)


def download_record_text(cohort_name, spec, attrition, selected, download_check,
                         sample_qc=None, out_dir=None, download_dir=None):
    """Render the cases-and-files record as text.

    The text is a pure function of its arguments: it holds no run timestamp, so
    two runs that reach the same result produce identical text. That is what
    makes :func:`write_if_changed` able to leave the file alone.

    Parameters
    ----------
    cohort_name : str
        Cohort label, e.g. ``"Lung_MALE"``.
    spec : gdc2maf.spec.FileSpec
        Which GDC files the cohort's MAFs are.
    attrition : pd.DataFrame
        Output of :func:`gdc2maf.reports.summarize_attrition`: one row per
        (step, reason).
    selected : pd.DataFrame
        The selected files, one per patient.
    download_check : pd.DataFrame
        Output of :func:`gdc2maf.download.check_downloaded_files`.
    sample_qc : pd.DataFrame or None, optional
        Per-sample QC from the merge, to report patients with and without
        variants. ``None`` leaves that section out.
    out_dir : str or None, optional
        Cohort output directory, read for the query and gdc-client records so
        the GDC data release and client version can be named.
    download_dir : str or None, optional
        Directory the files were downloaded into, named in the text.

    Returns
    -------
    str
        The record, ending in a newline.
    """
    query = _query_record(out_dir, cohort_name) if out_dir else {}
    client = _gdc_client_record(out_dir, cohort_name) if out_dir else {}

    lines = [
        "gdc2maf download record",
        "=======================",
        "",
        f"cohort:           {cohort_name}",
        f"projects:         {', '.join(query.get('projects', [])) or '-'}",
        f"sex at birth:     {', '.join(query.get('sex_at_birth') or []) or 'any'}",
        f"files:            {spec.title}",
    ]
    if query.get("gdc_data_release"):
        lines.append(
            f"GDC release:      {query['gdc_data_release']} "
            f"(queried {query.get('query_date', '-')})"
        )
    if client.get("version"):
        lines.append(
            f"gdc-client:       {client['version']} ({client.get('source', '-')})"
        )
    if download_dir:
        lines.append(f"downloaded into:  {download_dir}")

    lines += ["", "Cases", "-----"]
    for _, row in attrition.iterrows():
        if row["reason"] == "all cases":
            lines.append(f"  {int(row['n_remaining']):>7}  {row['step']}")
            continue
        if row["n_removed"] == 0:
            lines.append(f"  {'-':>7}  {row['step']}: nothing excluded")
            continue
        removed = f"-{int(row['n_removed'])}"
        lines.append(
            f"  {removed:>7}  {row['step']}: {row['reason']}"
            f"  (remaining {int(row['n_remaining'])})"
        )
    lines.append(f"  {int(attrition['n_remaining'].iloc[-1]):>7}  included")

    ok = download_check["status"] == "ok"
    lines += [
        "",
        "Files",
        "-----",
        f"  {len(selected):>7}  selected (one per patient)",
        f"  {int(ok.sum()):>7}  downloaded and verified (size + md5)",
        f"  {int((~ok).sum()):>7}  failed verification",
        f"  {_format_size(selected['file_size'].sum()):>7}  total size",
    ]
    if int((~ok).sum()):
        for status, n in download_check.loc[~ok, "status"].value_counts().items():
            lines.append(f"           {status}: {n}")

    if sample_qc is not None:
        with_variants = sample_qc["has_variants"].astype(bool)
        lines += [
            "",
            "Merged MAF",
            "----------",
            f"  {int(with_variants.sum()):>7}  patients with variants",
            f"  {int((~with_variants).sum()):>7}  patients with no variants",
        ]
        for reason, n in (
            sample_qc.loc[~with_variants, "no_variants_reason"].value_counts().items()
        ):
            lines.append(f"           {reason}: {n}")

    return "\n".join(lines) + "\n"


def write_if_changed(path, text, keep_superseded=True):
    """Write ``text`` to ``path`` only if the file does not already say it.

    Parameters
    ----------
    path : str
        File to write.
    text : str
        The content the file should have.
    keep_superseded : bool, optional
        When the file exists with different content, copy the old content to
        ``{stem}_superseded_{date}{suffix}`` before writing, so no earlier
        record is lost. A second change on the same day gets ``_2``, ``_3``…

    Returns
    -------
    dict
        ``status`` (``"unchanged"``, ``"created"`` or ``"updated"``), ``path``,
        and ``superseded`` (path of the kept copy, or ``None``).
    """
    if os.path.exists(path):
        with open(path) as f:
            if f.read() == text:
                logger.info(f"  record unchanged: {path}")
                return {"status": "unchanged", "path": path, "superseded": None}

    superseded = None
    if os.path.exists(path) and keep_superseded:
        stem, suffix = os.path.splitext(path)
        superseded = f"{stem}_superseded_{date.today().isoformat()}{suffix}"
        n = 2
        while os.path.exists(superseded):
            superseded = f"{stem}_superseded_{date.today().isoformat()}_{n}{suffix}"
            n += 1
        os.rename(path, superseded)
        logger.info(f"  previous record kept as {superseded}")

    status = "updated" if superseded else "created"
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as f:
        f.write(text)
    logger.info(f"  record {status}: {path}")
    return {"status": status, "path": path, "superseded": superseded}
