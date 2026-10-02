"""Locating a gdc-client, and what the platform table does and does not gate."""

from unittest import mock

import pytest

from gdc2maf.client import GDC_PLATFORM, ensure_gdc_client, install_gdc_client


@pytest.fixture
def fake_client(tmp_path):
    """An executable standing in for gdc-client, reporting version 2.3."""
    path = tmp_path / "gdc-client"
    path.write_text('#!/bin/sh\necho "2.3"\n')
    path.chmod(0o755)
    return str(path)


def test_the_gdc_publishes_builds_for_three_platforms():
    assert set(GDC_PLATFORM) == {"Linux", "Darwin", "Windows"}


def test_downloading_is_refused_where_no_build_is_published():
    with mock.patch("platform.system", return_value="FreeBSD"):
        with pytest.raises(ValueError) as exc:
            install_gdc_client()
    message = str(exc.value)
    assert "FreeBSD" in message
    # the error has to say what to do instead, not just refuse
    assert "--gdc-client" in message
    assert "gdc_client_path" in message


def test_a_supplied_client_is_used_on_a_platform_with_no_published_build(fake_client):
    # the platform table gates downloading a prebuilt release, not running
    with mock.patch("platform.system", return_value="FreeBSD"):
        info = ensure_gdc_client(gdc_client_path=fake_client)
    assert info["source"] == "user-provided"
    assert info["version"] == "2.3"
    assert info["path"] == fake_client


def test_a_supplied_client_records_its_own_md5(fake_client):
    from gdc2maf import file_md5

    info = ensure_gdc_client(gdc_client_path=fake_client)
    assert info["binary_md5"] == file_md5(fake_client)


def test_a_version_mismatch_warns_but_still_uses_the_client(fake_client, caplog):
    import logging

    with caplog.at_level(logging.WARNING, logger="gdc2maf.client"):
        info = ensure_gdc_client(version="2.0", gdc_client_path=fake_client)
    assert info["version"] == "2.3"
    assert "pinned version is 2.0" in caplog.text


def test_latest_does_not_warn_about_the_version(fake_client, caplog):
    import logging

    with caplog.at_level(logging.WARNING, logger="gdc2maf.client"):
        ensure_gdc_client(version="latest", gdc_client_path=fake_client)
    assert "pinned version" not in caplog.text


def test_a_missing_client_path_is_an_error(tmp_path):
    with pytest.raises(FileNotFoundError):
        ensure_gdc_client(gdc_client_path=str(tmp_path / "nope"))


def test_a_non_executable_client_path_is_an_error(tmp_path):
    path = tmp_path / "gdc-client"
    path.write_text("not executable")
    path.chmod(0o644)
    with pytest.raises(FileNotFoundError):
        ensure_gdc_client(gdc_client_path=str(path))
