"""TCGA aliquot barcode parsing."""

import pandas as pd

from gdc2maf.tcga import parse_tcga_aliquot_barcode


def test_every_component_is_split_out():
    parsed = parse_tcga_aliquot_barcode(pd.Series(["TCGA-44-6147-01A-11D-A271-08"]))
    row = parsed.iloc[0]
    assert row["patient"] == "TCGA-44-6147"
    assert row["sample_type_code"] == 1
    assert row["vial"] == "A"
    assert row["portion"] == 11
    assert row["analyte"] == "D"
    assert row["plate"] == "A271"
    assert row["center"] == "08"


def test_wga_analyte_and_normal_code_are_read():
    parsed = parse_tcga_aliquot_barcode(
        pd.Series(["TCGA-44-6147-10A-01W-A271-08"])
    )
    assert parsed.iloc[0]["analyte"] == "W"
    assert parsed.iloc[0]["sample_type_code"] == 10


def test_index_is_preserved():
    barcodes = pd.Series(["TCGA-44-6147-01A-11D-A271-08"], index=[7])
    assert list(parse_tcga_aliquot_barcode(barcodes).index) == [7]
