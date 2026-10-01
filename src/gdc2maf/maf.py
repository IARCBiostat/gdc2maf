"""Merge a cohort's downloaded MAFs into one MAF, with QC.

The merge is strict about anything that would silently corrupt the result (a
file that failed its download check, a header that does not match, a row on the
wrong reference build, a barcode that is not the selected aliquot, a duplicated
variant) and raises rather than guessing.

It is deliberately *not* strict about biology: samples with an implausible
mutation count, excess low-VAF calls or weak caller support are flagged with
the evidence in the per-sample QC table and left in the MAF, because whether to
drop them is the analyst's decision, not the loader's.
"""

import gzip
import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

EXPECTED_BUILD = "GRCh38"
WGA_ANALYTES = ["W", "X"]
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
SBS_CHANGES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
# Fixed report order; values not listed here are added after these
VARIANT_TYPES = ["SNP", "DNP", "TNP", "ONP", "INS", "DEL"]
GDC_FILTER_FLAGS = ["NonExonic", "gdc_pon", "common_in_gnomAD"]
NO_VARIANTS_EMPTY = "empty GDC MAF"
NO_VARIANTS_FILTERED = "only variants with removed GDC_FILTER flags"

# Artefact indicators (per sample). A sample is flagged for review when its
# count is extreme, or it is in the top 5% of counts and shows signs of
# technical noise.
LOW_VAF = 0.1
MAX_FRAC_LOW_VAF = 0.25  # >= 25% of variants with VAF < 0.1
MAX_FRAC_WEAK_CALLER_SUPPORT = 0.35  # >= 35% of variants called by <= 2 callers
HIGH_COUNT_IQR_FACTOR = 3  # log10(count) > Q3 + 3 * IQR within the cohort


def read_gdc_maf(path):
    """Read one GDC MAF (gzipped), skipping ``#`` comment lines.

    All values are kept as text exactly as in the file (no type conversion,
    empty strings stay empty).

    Returns
    -------
    comments : list of str
        The ``#`` lines.
    header : list of str
        Column names.
    maf : pd.DataFrame
        Variant rows (may be empty).
    """
    with gzip.open(path, "rt") as f:
        lines = f.read().splitlines()
    comments = [line for line in lines if line.startswith("#")]
    body = [line for line in lines if not line.startswith("#")]
    header = body[0].split("\t")
    rows = [
        line.split("\t")
        for line in body[1:]
        if line and not line.startswith("Hugo_Symbol\t")
    ]
    return comments, header, pd.DataFrame(rows, columns=header, dtype=str)


def sample_qc_metrics(maf):
    """Per-sample metrics used to spot extreme or artefactual mutation counts.

    Parameters
    ----------
    maf : pd.DataFrame
        Merged MAF (text values).

    Returns
    -------
    pd.DataFrame
        Indexed by ``Tumor_Sample_Barcode``: ``n_variants``, ``n_snv``,
        ``n_indel``, ``median_vaf``, ``frac_vaf_lt_0.1``,
        ``frac_weak_caller_support`` (<= 2 callers), ``n_gdc_filter_flagged``
        (variants with any ``GDC_FILTER`` flag, kept), ``frac_nonexonic`` and
        the fraction of SNVs in each of the six pyrimidine-based
        substitution classes.
    """
    vaf = pd.to_numeric(maf["t_alt_count"]) / pd.to_numeric(maf["t_depth"])
    df = pd.DataFrame(
        {
            "Tumor_Sample_Barcode": maf["Tumor_Sample_Barcode"],
            "is_snv": maf["Variant_Type"] == "SNP",
            "is_indel": maf["Variant_Type"].isin(["INS", "DEL"]),
            "vaf": vaf,
            "low_vaf": vaf < LOW_VAF,
            "weak_callers": maf["callers"].str.count(";") <= 1,
            "nonexonic": maf["GDC_FILTER"].str.contains("NonExonic"),
            "gdc_filter_flagged": maf["GDC_FILTER"] != "",
        }
    )
    g = df.groupby("Tumor_Sample_Barcode")
    metrics = pd.DataFrame(
        {
            "n_variants": g.size(),
            "n_snv": g["is_snv"].sum(),
            "n_indel": g["is_indel"].sum(),
            "median_vaf": g["vaf"].median(),
            "frac_vaf_lt_0.1": g["low_vaf"].mean(),
            "frac_weak_caller_support": g["weak_callers"].mean(),
            "n_gdc_filter_flagged": g["gdc_filter_flagged"].sum(),
            "frac_nonexonic": g["nonexonic"].mean(),
        }
    )

    snv = maf[df["is_snv"]]
    ref, alt = snv["Reference_Allele"], snv["Tumor_Seq_Allele2"]
    purine = ref.isin(["A", "G"])
    change = np.where(
        purine, ref.map(COMPLEMENT) + ">" + alt.map(COMPLEMENT), ref + ">" + alt
    )
    spectrum = (
        pd.crosstab(snv["Tumor_Sample_Barcode"], change, normalize="index")
        .reindex(columns=SBS_CHANGES, fill_value=0)
        .add_prefix("frac_")
    )
    return metrics.join(spectrum).fillna({f"frac_{c}": 0 for c in SBS_CHANGES})


def flag_samples(sample_qc):
    """Flag samples with an extreme count and/or artefact indicators.

    - ``high_count``: log10(n_variants) > Q3 + 3 * IQR of the cohort.
    - ``low_vaf_excess``: >= 25% of variants with VAF < 0.1.
    - ``weak_caller_support``: >= 35% of variants called by <= 2 callers.
    - ``wga``: tumour DNA is whole-genome amplified.

    ``review`` is True for high-count samples and for samples in the top 5%
    of counts with low-VAF excess or weak caller support. WGA is reported but
    does not trigger review on its own, since whole-genome amplification is a
    property of the library rather than evidence of a bad call; handle those
    patients with a sensitivity analysis instead.

    Parameters
    ----------
    sample_qc : pd.DataFrame
        Output of :func:`sample_qc_metrics`, with ``wga`` added.

    Returns
    -------
    pd.DataFrame
        ``sample_qc`` with ``high_count_threshold``, ``high_count``,
        ``low_vaf_excess``, ``weak_caller_support`` and ``review``.
    """
    qc = sample_qc.copy()
    logn = np.log10(qc["n_variants"])
    q1, q3 = logn.quantile([0.25, 0.75])
    qc["high_count_threshold"] = int(
        round(10 ** (q3 + HIGH_COUNT_IQR_FACTOR * (q3 - q1)))
    )
    qc["high_count"] = logn > q3 + HIGH_COUNT_IQR_FACTOR * (q3 - q1)
    qc["low_vaf_excess"] = qc["frac_vaf_lt_0.1"] >= MAX_FRAC_LOW_VAF
    qc["weak_caller_support"] = (
        qc["frac_weak_caller_support"] >= MAX_FRAC_WEAK_CALLER_SUPPORT
    )
    indicators = qc[["low_vaf_excess", "weak_caller_support"]].any(axis=1)
    top5 = qc["n_variants"] >= qc["n_variants"].quantile(0.95)
    qc["review"] = qc["high_count"] | (top5 & indicators)
    return qc


def merge_mafs(
    cohort_name,
    selected,
    download_check,
    out_dir,
    remove_gdc_filter_flags=(),
    maf_name=None,
):
    """Merge a cohort's selected MAFs into one MAF and run QC.

    Hard checks (raise ``ValueError``): every selected file passed the
    download check; all files have the same header; every row is
    ``NCBI_Build == GRCh38``; each file's ``Tumor_Sample_Barcode`` and
    ``Matched_Norm_Sample_Barcode`` are the selected aliquots; no duplicate
    variant within a sample.

    Reported (not fatal): files with no variants; ``Variant_Type`` counts;
    ``GDC_FILTER`` counts per flag; per-sample metrics and flags for extreme
    counts. By default all ``GDC_FILTER``-flagged variants are kept; flags
    listed in ``remove_gdc_filter_flags`` are removed.

    Two columns are appended to the MAF: ``tumor_analyte`` and
    ``tumor_is_wga`` (from the tumour aliquot barcode).

    Parameters
    ----------
    cohort_name : str
        Cohort label, e.g. ``"Lung_MALE"``.
    selected : pd.DataFrame
        Selected files from
        :func:`gdc2maf.selection.select_one_file_per_case`.
    download_check : pd.DataFrame
        Output of :func:`gdc2maf.download.download_files`.
    out_dir : str
        Directory for the merged MAF and the QC tables.
    remove_gdc_filter_flags : iterable of str, optional
        ``GDC_FILTER`` flags (e.g. ``"NonExonic"``, ``"gdc_pon"``,
        ``"common_in_gnomAD"``) whose variants are removed before the MAF is
        written and QC metrics are computed. A variant is removed if any of
        its flags is listed. Empty (default) keeps all variants. Patients left
        with no variants get ``no_variants_reason`` =
        ``"only variants with removed GDC_FILTER flags"``.
    maf_name : str or None, optional
        File name for the merged MAF inside ``out_dir``. ``None`` uses
        ``{cohort_name}.maf``.

    Returns
    -------
    maf : pd.DataFrame
        The merged MAF.
    sample_qc : pd.DataFrame
        One row per selected patient (including those with no variants):
        metrics, flags, ``has_variants`` and ``no_variants_reason``.
    qc_summary : dict
        Cohort-level QC counts for the cross-cohort report.
    """
    not_ok = download_check[download_check["status"] != "ok"]
    if len(not_ok):
        raise ValueError(
            f"{cohort_name}: {len(not_ok)} selected files failed the download check"
        )
    paths = download_check.set_index("file_id")["path"]

    header, frames, n_comment_lines, empty = None, [], set(), []
    for _, f in selected.iterrows():
        comments, file_header, df = read_gdc_maf(paths[f["file_id"]])
        n_comment_lines.add(len(comments))
        if header is None:
            header = file_header
        elif file_header != header:
            raise ValueError(
                f"{cohort_name}: header of {f['file_id']} differs from the first file"
            )
        if df.empty:
            empty.append(f["submitter_id"])
            continue
        if set(df["Tumor_Sample_Barcode"]) != {f["tumor_aliquot_barcode"]}:
            raise ValueError(
                f"{cohort_name}: {f['file_id']} tumour barcode "
                f"is not the selected aliquot"
            )
        if set(df["Matched_Norm_Sample_Barcode"]) != {f["normal_aliquot_barcode"]}:
            raise ValueError(
                f"{cohort_name}: {f['file_id']} normal barcode "
                f"is not the selected aliquot"
            )
        frames.append(df)
    maf = pd.concat(frames, ignore_index=True)

    builds = set(maf["NCBI_Build"])
    if builds != {EXPECTED_BUILD}:
        raise ValueError(
            f"{cohort_name}: NCBI_Build values {builds}, expected {EXPECTED_BUILD}"
        )
    variant_key = [
        "Tumor_Sample_Barcode",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
    ]
    n_duplicates = int(maf.duplicated(variant_key).sum())
    if n_duplicates:
        raise ValueError(
            f"{cohort_name}: {n_duplicates} duplicate variants within samples"
        )

    # GDC_FILTER can hold several flags separated by ";"; count variants per flag
    n_variants_gdc = len(maf)
    flag_lists = maf["GDC_FILTER"].str.split(";")
    flag_counts = flag_lists[maf["GDC_FILTER"] != ""].explode().value_counts()
    flag_counts = flag_counts.reindex(
        GDC_FILTER_FLAGS + sorted(set(flag_counts.index) - set(GDC_FILTER_FLAGS)),
        fill_value=0,
    )
    n_flagged_gdc = int((maf["GDC_FILTER"] != "").sum())

    remove_gdc_filter_flags = list(dict.fromkeys(remove_gdc_filter_flags))
    samples_before = set(maf["Tumor_Sample_Barcode"])
    to_remove = flag_lists.apply(
        lambda flags: any(f in remove_gdc_filter_flags for f in flags)
    )
    # removed variants per listed flag (a variant with several listed flags
    # counts for each)
    removed_per_flag = {
        flag: int(
            flag_lists[to_remove].apply(lambda flags, flag=flag: flag in flags).sum()
        )
        for flag in remove_gdc_filter_flags
    }
    n_removed = int(to_remove.sum())
    maf = maf[~to_remove].reset_index(drop=True)
    only_removed_flags = samples_before - set(maf["Tumor_Sample_Barcode"])

    analyte = maf["Tumor_Sample_Barcode"].str.split("-").str[4].str[-1]
    maf["tumor_analyte"] = analyte
    maf["tumor_is_wga"] = analyte.isin(WGA_ANALYTES)
    maf_path = f"{out_dir}/{maf_name or f'{cohort_name}.maf'}"
    maf.to_csv(maf_path, sep="\t", index=False)

    # Per-patient QC, including patients whose MAF has no variants
    patients = selected.set_index("tumor_aliquot_barcode")[
        ["submitter_id", "case_id", "tumor_analyte"]
    ]
    sample_qc = patients.join(sample_qc_metrics(maf), how="left")
    sample_qc.index.name = "Tumor_Sample_Barcode"
    sample_qc["has_variants"] = sample_qc["n_variants"].notna()
    sample_qc["no_variants_reason"] = np.where(
        sample_qc["has_variants"],
        "",
        np.where(
            sample_qc.index.isin(only_removed_flags),
            NO_VARIANTS_FILTERED,
            NO_VARIANTS_EMPTY,
        ),
    )
    sample_qc["n_variants"] = sample_qc["n_variants"].fillna(0).astype(int)
    sample_qc["wga"] = sample_qc["tumor_analyte"].isin(WGA_ANALYTES)
    flagged = flag_samples(sample_qc[sample_qc["has_variants"]])
    sample_qc = sample_qc.join(
        flagged[
            [
                "high_count_threshold",
                "high_count",
                "low_vaf_excess",
                "weak_caller_support",
                "review",
            ]
        ]
    )
    sample_qc = sample_qc.sort_values("n_variants", ascending=False).reset_index()
    sample_qc.to_csv(f"{out_dir}/{cohort_name}_sample_qc.tsv", sep="\t", index=False)

    n_samples = maf["Tumor_Sample_Barcode"].nunique()
    n_with_variants = len(selected) - len(empty) - len(only_removed_flags)
    variant_types = maf["Variant_Type"].value_counts()
    variant_types = variant_types.reindex(
        VARIANT_TYPES + sorted(set(variant_types.index) - set(VARIANT_TYPES)),
        fill_value=0,
    )
    flagged_rows = maf["GDC_FILTER"] != ""
    pct_flagged = round(100 * flagged_rows.mean(), 2)

    qc_summary = {
        "Cohort": cohort_name,
        "Selected patients (files)": len(selected),
        "Comment lines per file": "/".join(str(n) for n in sorted(n_comment_lines)),
        "Files with no variants": len(empty),
        "NCBI_Build": ";".join(sorted(builds)),
        "Duplicate variants": n_duplicates,
        "Variants in GDC MAFs": n_variants_gdc,
        "GDC_FILTER flagged in GDC MAFs": n_flagged_gdc,
        **{f"  GDC_FILTER {k}": int(v) for k, v in flag_counts.items()},
        "GDC_FILTER flags removed": ";".join(remove_gdc_filter_flags)
        or "none (all kept)",
        "Variants removed (GDC_FILTER)": n_removed,
        **{f"  removed: {k}": v for k, v in removed_per_flag.items()},
        "Patients left with no variants after removal": len(only_removed_flags),
        "Variants in merged MAF": len(maf),
        **{f"  Variant_Type {k}": int(v) for k, v in variant_types.items()},
        "GDC_FILTER flagged, kept in merged MAF": int(flagged_rows.sum()),
        "GDC_FILTER flagged, kept (%)": pct_flagged,
        "Unique Tumor_Sample_Barcode": n_samples,
        "Barcodes = selected patients with variants": "yes"
        if n_samples == n_with_variants
        else "NO",
        "Median variants per sample": float(
            sample_qc.loc[sample_qc["has_variants"], "n_variants"].median()
        ),
        "Max variants per sample": int(sample_qc["n_variants"].max()),
        "High-count threshold": int(flagged["high_count_threshold"].iloc[0]),
        "High-count samples": int(flagged["high_count"].sum()),
        "Samples flagged for review": int(flagged["review"].sum()),
    }
    logger.info(
        f"{cohort_name}: merged {len(selected) - len(empty)} MAFs "
        f"({len(maf)} variants, {n_samples} samples; {len(empty)} files with "
        f"no variants: {', '.join(empty) or '-'}) -> {maf_path}"
    )
    logger.info(
        f"  GDC_FILTER-flagged variants in GDC MAFs: {n_flagged_gdc} "
        f"({round(100 * n_flagged_gdc / n_variants_gdc, 2)}%): "
        + ", ".join(f"{k} {int(v)}" for k, v in flag_counts.items())
    )
    if remove_gdc_filter_flags:
        logger.info(
            f"  GDC_FILTER flags removed: {n_removed} variants ("
            + ", ".join(f"{k} {v}" for k, v in removed_per_flag.items())
            + f"); patients left with no variants: {len(only_removed_flags)}"
        )
    else:
        logger.info("  GDC_FILTER flags removed: none (all flagged variants kept)")
    logger.info(
        f"  GDC_FILTER-flagged variants kept in merged MAF: "
        f"{int(flagged_rows.sum())} ({pct_flagged}%)"
    )
    logger.info(
        f"  samples flagged for review: {int(flagged['review'].sum())} "
        f"(see samples_for_review.tsv)"
    )
    return maf, sample_qc, qc_summary
