"""TCGA aliquot barcode parsing.

The barcode carries the facts the file-selection rule needs and that the GDC
file metadata does not expose directly: the sample type and vial, the portion,
the analyte (``D`` native DNA, ``W``/``X`` whole-genome amplified), the
sequencing plate and the centre.
"""

import pandas as pd


def parse_tcga_aliquot_barcode(barcodes):
    """Split TCGA aliquot barcodes into their components.

    ``TCGA-44-6147-01A-11D-A271-08`` is
    project-TSS-participant-{sample type}{vial}-{portion}{analyte}-plate-center.

    Parameters
    ----------
    barcodes : pd.Series of str
        Aliquot barcodes.

    Returns
    -------
    pd.DataFrame
        Same index as ``barcodes`` with columns ``patient``,
        ``sample_type_code`` (int; 01 primary, 06 metastatic, 10 blood normal,
        11 solid-tissue normal...), ``vial``, ``portion`` (int), ``analyte``
        (``D`` native DNA, ``W``/``X`` whole-genome amplified), ``plate`` and
        ``center``.
    """
    parts = barcodes.str.split("-", expand=True)
    return pd.DataFrame(
        {
            "patient": parts[0] + "-" + parts[1] + "-" + parts[2],
            "sample_type_code": parts[3].str[:2].astype(int),
            "vial": parts[3].str[2],
            "portion": parts[4].str[:2].astype(int),
            "analyte": parts[4].str[2],
            "plate": parts[5],
            "center": parts[6],
        },
        index=barcodes.index,
    )
