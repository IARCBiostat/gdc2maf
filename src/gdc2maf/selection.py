"""Choose one MAF file per patient.

A patient with several MAF files would otherwise contribute several times to
the merged MAF. The rule here ("barcode rule") sorts a patient's files by
:data:`SELECTION_KEYS` and keeps the first, recording for every dropped file
*which* key decided against it, so the choice is auditable rather than
arbitrary.

The keys read TCGA aliquot barcodes (see :mod:`gdc2maf.tcga`), so this rule
applies to TCGA projects. For other projects, sort ``maf_files`` yourself and
keep one row per ``case_id``.
"""

import logging

from .tcga import parse_tcga_aliquot_barcode

logger = logging.getLogger(__name__)

# Primary-tumour sample type codes: 01 solid primary, 03 blood-derived primary,
# 09 bone-marrow primary.
PRIMARY_TUMOUR_CODES = [1, 3, 9]
NATIVE_DNA_ANALYTES = ["D"]
# 10 blood-derived normal before 11 solid-tissue normal: adjacent normal tissue
# can carry tumour cells, which hides true somatic variants.
NORMAL_CODE_RANK = {10: 0, 11: 1}

# Ordered selection keys: (column, ascending, label). Files of a patient are
# sorted by these keys and the first file is kept. The last two keys make the
# order total, so the choice never depends on input row order.
SELECTION_KEYS = [
    ("rank_primary", True, "non-primary tumour (primary available)"),
    ("rank_ffpe", True, "FFPE tumour (non-FFPE available)"),
    ("rank_native", True, "WGA tumour DNA (native DNA available)"),
    ("tumor_vial", True, "later tumour vial"),
    ("tumor_plate", False, "earlier sequencing plate"),
    ("tumor_portion", True, "higher tumour portion"),
    ("rank_normal", True, "solid-tissue normal (blood normal available)"),
    ("tumor_aliquot_barcode", True, "tumour aliquot barcode order"),
    ("normal_aliquot_barcode", True, "normal aliquot barcode order"),
]


def add_selection_ranks(maf_files):
    """Add the barcode fields and rank columns used by :data:`SELECTION_KEYS`.

    Parameters
    ----------
    maf_files : pd.DataFrame
        Output of :func:`gdc2maf.files.list_maf_files`.

    Returns
    -------
    pd.DataFrame
        Copy of ``maf_files`` with ``tumor_*`` / ``normal_*`` barcode fields
        and ``rank_*`` columns (0 = preferred).
    """
    files = maf_files.copy()
    tumor = parse_tcga_aliquot_barcode(files["tumor_aliquot_barcode"]).add_prefix(
        "tumor_"
    )
    normal = parse_tcga_aliquot_barcode(files["normal_aliquot_barcode"]).add_prefix(
        "normal_"
    )
    files = files.join(tumor.drop(columns="tumor_patient")).join(
        normal.drop(columns="normal_patient")
    )

    files["rank_primary"] = (
        ~files["tumor_sample_type_code"].isin(PRIMARY_TUMOUR_CODES)
    ).astype(int)
    # FFPE samples carry sample type code 01, so use the GDC sample type label
    files["rank_ffpe"] = files["tumor_sample_type"].str.contains("FFPE").astype(int)
    files["rank_native"] = (~files["tumor_analyte"].isin(NATIVE_DNA_ANALYTES)).astype(
        int
    )
    files["rank_normal"] = (
        files["normal_sample_type_code"].map(NORMAL_CODE_RANK).fillna(2).astype(int)
    )
    return files


def drop_reason(dropped, selected):
    """Name the first selection key on which ``dropped`` lost to ``selected``."""
    for col, _, label in SELECTION_KEYS:
        if dropped[col] != selected[col]:
            return label
    raise ValueError(
        f"{dropped['file_id']}: identical selection keys to the selected file"
    )


def select_one_file_per_case(cohort_name, maf_files, out_dir):
    """Select one MAF file per patient with the barcode rule (Option A).

    Files are sorted by :data:`SELECTION_KEYS` within each patient and the
    first is kept:

    1. primary tumour over metastatic / recurrent / new primary
    2. non-FFPE over FFPE tumour
    3. native DNA (``D``) over whole-genome amplified (``W``/``X``)
    4. lower vial (``A`` before ``B``)
    5. later plate (plate codes compared as text, e.g. ``A27T`` > ``A271`` > ``1753``)
    6. lower portion
    7. blood-derived normal over solid-tissue normal
    8. tumour, then normal aliquot barcode, as a final deterministic tie-break

    Parameters
    ----------
    cohort_name : str
        Cohort label used in file names, e.g. ``"Lung_MALE"``.
    maf_files : pd.DataFrame
        Output of :func:`gdc2maf.files.list_maf_files`.
    out_dir : str
        Directory for ``{cohort_name}_files_selected.tsv`` and
        ``{cohort_name}_files_dropped.tsv``.

    Returns
    -------
    selected : pd.DataFrame
        One row per patient: the chosen file with barcode fields, ranks and
        ``n_files`` for the patient.
    dropped : pd.DataFrame
        Every other file, with ``selected_file_id`` and ``drop_reason``.
    """
    files = add_selection_ranks(maf_files)
    files["n_files"] = files.groupby("case_id")["file_id"].transform("size")

    cols = ["case_id"] + [col for col, _, _ in SELECTION_KEYS]
    ascending = [True] + [asc for _, asc, _ in SELECTION_KEYS]
    files = files.sort_values(cols, ascending=ascending, kind="stable")

    is_selected = ~files.duplicated("case_id", keep="first")
    selected = files[is_selected].reset_index(drop=True)
    dropped = files[~is_selected].copy()

    selected_by_case = selected.set_index("case_id")
    dropped["selected_file_id"] = dropped["case_id"].map(selected_by_case["file_id"])
    dropped["drop_reason"] = [
        drop_reason(row, selected_by_case.loc[row["case_id"]])
        for _, row in dropped.iterrows()
    ]
    dropped = dropped.reset_index(drop=True)

    selected.to_csv(
        f"{out_dir}/{cohort_name}_files_selected.tsv", sep="\t", index=False
    )
    dropped.to_csv(f"{out_dir}/{cohort_name}_files_dropped.tsv", sep="\t", index=False)
    logger.info(
        f"{cohort_name}: selected {len(selected)} files, dropped {len(dropped)} "
        f"-> {out_dir}/{cohort_name}_files_selected.tsv"
    )
    return selected, dropped
