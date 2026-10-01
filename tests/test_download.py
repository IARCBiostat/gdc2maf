"""Manifest writing and the size + md5 verification of downloaded files."""

import gzip
import hashlib

import pandas as pd

from gdc2maf.download import check_downloaded_files, file_md5, write_manifest


def test_manifest_has_the_gdc_column_names(tmp_path, maf_files):
    path = tmp_path / "manifest.txt"
    write_manifest(maf_files, str(path))
    manifest = pd.read_csv(path, sep="\t")
    assert list(manifest.columns) == ["id", "filename", "md5", "size", "state"]
    assert len(manifest) == len(maf_files)


def test_file_md5_matches_hashlib(tmp_path):
    path = tmp_path / "f.bin"
    payload = b"some bytes" * 1000
    path.write_bytes(payload)
    assert file_md5(str(path)) == hashlib.md5(payload).hexdigest()


def _one_file(tmp_path, content=b"maf content"):
    """Write one file where gdc-client would put it; return its metadata row."""
    file_dir = tmp_path / "fid"
    file_dir.mkdir()
    path = file_dir / "f.maf.gz"
    with gzip.open(path, "wb") as f:
        f.write(content)
    return pd.DataFrame(
        [
            {
                "case_id": "case-1",
                "submitter_id": "TCGA-AA-0001",
                "file_id": "fid",
                "file_name": "f.maf.gz",
                "file_size": path.stat().st_size,
                "md5sum": file_md5(str(path)),
            }
        ]
    ), path


def test_a_valid_file_checks_out_ok(tmp_path):
    files, _ = _one_file(tmp_path)
    check = check_downloaded_files(files, str(tmp_path))
    assert list(check["status"]) == ["ok"]


def test_a_missing_file_is_reported(tmp_path):
    files, path = _one_file(tmp_path)
    path.unlink()
    check = check_downloaded_files(files, str(tmp_path))
    assert check.loc[0, "status"] == "missing"


def test_a_corrupt_file_is_reported(tmp_path):
    files, path = _one_file(tmp_path)
    path.write_bytes(b"truncated")
    check = check_downloaded_files(files, str(tmp_path))
    assert check.loc[0, "status"] in {"size mismatch", "md5 mismatch"}


def test_check_reports_the_path_it_looked_at(tmp_path):
    files, path = _one_file(tmp_path)
    check = check_downloaded_files(files, str(tmp_path))
    assert check.loc[0, "path"] == str(path)
