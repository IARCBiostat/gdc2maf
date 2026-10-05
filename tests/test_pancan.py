"""The PanCanAtlas supplements: pinned downloads and the Do_not_use flag."""

import hashlib

import pytest

from gdc2maf.pancan import (
    DO_NOT_USE_COLUMN,
    PANCAN_FILES,
    do_not_use_patients,
    fetch_pancan_file,
    load_quality_annotations,
    write_pancan_record,
)

QUALITY_TSV = (
    "patient_barcode\taliquot_barcode\tcancer type\tplatform\tDo_not_use\n"
    "TCGA-AA-0001\tTCGA-AA-0001-01A-11D-A271-08\tLUAD\tWXS\tFalse\n"
    "TCGA-AA-0002\tTCGA-AA-0002-01A-11D-A271-08\tLUAD\tWXS\tTrue\n"
    "TCGA-AA-0003\tTCGA-AA-0003-01A-11D-A271-08\tLUAD\tSNP6\tTrue\n"
    "TCGA-AA-0003\tTCGA-AA-0003-01A-11D-A271-08\tLUAD\tWXS\tFalse\n"
)


@pytest.fixture
def quality_file(tmp_path):
    """A small stand-in for merged_sample_quality_annotations.tsv."""
    path = tmp_path / "merged_sample_quality_annotations.tsv"
    path.write_text(QUALITY_TSV)
    return path


def test_the_registry_pins_a_uuid_and_an_md5_for_each_file():
    assert set(PANCAN_FILES) == {"clinical", "quality_annotations"}
    for name, (file_name, uuid, md5, description) in PANCAN_FILES.items():
        assert file_name.endswith(".tsv"), name
        assert len(uuid) == 36, name
        assert len(md5) == 32, name
        assert description, name


def test_an_unknown_file_name_is_an_error(tmp_path):
    with pytest.raises(KeyError, match="unknown PanCanAtlas file"):
        fetch_pancan_file("nonsense", str(tmp_path))


def test_an_existing_file_with_the_right_md5_is_not_downloaded(tmp_path, monkeypatch):
    payload = b"patient_barcode\tDo_not_use\nTCGA-AA-0001\tFalse\n"
    name, (file_name, uuid, _, description) = "clinical", PANCAN_FILES["clinical"]
    path = tmp_path / file_name
    path.write_bytes(payload)
    # pin the md5 to what is on disk, so no download is needed
    monkeypatch.setitem(
        PANCAN_FILES, name,
        (file_name, uuid, hashlib.md5(payload).hexdigest(), description),
    )

    def fail(*a, **k):
        raise AssertionError("should not download a valid existing file")

    monkeypatch.setattr("gdc2maf.pancan._download", fail)
    record = fetch_pancan_file(name, str(tmp_path))
    assert record["source"] == "existing"
    assert record["md5"] == record["expected_md5"]


def test_a_changed_upstream_file_fails_loudly(tmp_path, monkeypatch):
    file_name = PANCAN_FILES["clinical"][0]

    def serve_something_else(url, dest):
        with open(dest, "wb") as f:
            f.write(b"not the pinned file")

    monkeypatch.setattr("gdc2maf.pancan._download", serve_something_else)
    with pytest.raises(RuntimeError, match="md5 mismatch"):
        fetch_pancan_file("clinical", str(tmp_path))
    assert (tmp_path / file_name).exists()  # kept, so it can be inspected


def test_the_md5_check_can_be_waived(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "gdc2maf.pancan._download",
        lambda url, dest: open(dest, "wb").write(b"whatever"),
    )
    record = fetch_pancan_file("clinical", str(tmp_path), check_md5=False)
    assert record["md5"] != record["expected_md5"]


def test_the_record_is_json_serialisable(tmp_path, quality_file, monkeypatch):
    monkeypatch.setitem(
        PANCAN_FILES, "quality_annotations",
        (quality_file.name, "u" * 36,
         hashlib.md5(QUALITY_TSV.encode()).hexdigest(), "d"),
    )
    record = fetch_pancan_file("quality_annotations", str(quality_file.parent))
    out = write_pancan_record(record, str(tmp_path / "record.json"))
    import json

    assert json.loads(open(out).read())["source"] == "existing"


def test_do_not_use_is_read_as_a_boolean(quality_file):
    annotations = load_quality_annotations(str(quality_file))
    assert annotations[DO_NOT_USE_COLUMN].dtype == bool
    assert annotations[DO_NOT_USE_COLUMN].sum() == 2


def test_a_patient_flagged_on_any_platform_is_flagged(quality_file):
    # TCGA-AA-0003 is flagged on SNP6 but not on WXS; the conservative reading
    # the TCGA working groups applied flags the patient
    flagged = do_not_use_patients(load_quality_annotations(str(quality_file)))
    assert flagged == {"TCGA-AA-0002", "TCGA-AA-0003"}


def test_flagging_can_be_narrowed_to_one_platform(quality_file):
    annotations = load_quality_annotations(str(quality_file))
    assert do_not_use_patients(annotations, platform="WXS") == {"TCGA-AA-0002"}
    assert do_not_use_patients(annotations, platform="SNP6") == {"TCGA-AA-0003"}


def test_the_download_is_streamed_with_requests(tmp_path, monkeypatch):
    """The GDC closes the connection on urllib's default user agent."""
    import gdc2maf.pancan as pancan

    class FakeResponse:
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def raise_for_status(self):
            pass

        def iter_content(self, chunk_size):
            yield b"first "
            yield b"second"

    calls = {}

    def fake_get(url, **kwargs):
        calls.update(url=url, **kwargs)
        return FakeResponse()

    monkeypatch.setattr(pancan.requests, "get", fake_get)
    path = tmp_path / "out.tsv"
    pancan._download("https://example/data/x", str(path))
    assert path.read_bytes() == b"first second"
    assert calls["stream"] is True
    assert calls["timeout"] == 300
