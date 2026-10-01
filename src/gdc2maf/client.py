"""Find or install the GDC Data Transfer Tool (gdc-client).

The client is kept in ``{install_dir}/{version}/gdc-client``. If it is not
there, the release is downloaded from the GDC and md5-checked against
:data:`PINNED_RELEASES`, so a run never silently uses an unverified binary.
Pin a version for reproducible results; ``"latest"`` resolves the newest
release from the GDC download page.
"""

import json
import logging
import os
import platform
import re
import stat
import subprocess
import urllib.request
import zipfile
from datetime import date

from .download import file_md5

logger = logging.getLogger(__name__)

GDC_TOOL_PAGE = "https://gdc.cancer.gov/access-data/gdc-data-transfer-tool"
GDC_FILE_URL = "https://gdc.cancer.gov/system/files/public/file"

# Pinned releases: the version used for the published results.
# (zip name, md5 of the zip as published on GDC_TOOL_PAGE, md5 of the executable
# inside it; None where not checked).
PINNED_RELEASES = {
    "2.3": {
        "Linux": (
            "gdc-client_2.3_Ubuntu_x64-py3.8-ubuntu-20.04.zip",
            "18591d74de07cdcd396dab71c52663da",
            "45885253a71abdac17ca878abb413c30",
        ),
        "Darwin": (
            "gdc-client_2.3_OSX_x64-py3.8-macos-14.zip",
            "56cca3594fa5fb47bc8297f5b6fd0e20",
            None,
        ),
        "Windows": (
            "gdc-client_2.3_Windows_x64-py3.8-windows-2019.zip",
            "525ce44bb5f3f0624066b906c7dbdaf4",
            None,
        ),
    },
}
# Platform label used in GDC zip names
GDC_PLATFORM = {"Linux": "Ubuntu", "Darwin": "OSX", "Windows": "Windows"}


def gdc_client_version(gdc_client):
    """Return the version string printed by ``gdc-client --version``."""
    result = subprocess.run(
        [gdc_client, "--version"], capture_output=True, text=True, check=True
    )
    return result.stdout.strip()


def latest_gdc_client_release(system):
    """Find the newest gdc-client zip for ``system`` on the GDC tool page.

    Parameters
    ----------
    system : str
        ``platform.system()`` value: ``"Linux"``, ``"Darwin"`` or ``"Windows"``.

    Returns
    -------
    version : str
    zip_name : str
    md5 : str or None
        md5 of the zip as listed on the page; ``None`` if not listed.
    """
    html = (
        urllib.request.urlopen(GDC_TOOL_PAGE, timeout=60)
        .read()
        .decode("utf-8", "replace")
    )
    label = GDC_PLATFORM[system]
    zips = set(
        re.findall(rf"(gdc-client_([0-9][0-9.]*)_{label}_x64[^\"'<>\s]*\.zip)", html)
    )
    if not zips:
        raise RuntimeError(
            f"No gdc-client release for {system} found on {GDC_TOOL_PAGE}"
        )
    # For macOS several builds exist; prefer the newest macOS build
    zip_name, version = max(
        zips, key=lambda z: ([int(x) for x in z[1].split(".")], z[0])
    )

    text = re.sub(r"<[^>]+>", " ", html)
    md5 = None
    for name, checksum in re.findall(
        r"(gdc-client_[\w.\-]+)\s+md5sum:\s*([0-9a-f]{32})", text
    ):
        # page lists "gdc-client_2.3.0_Ubuntu_x64-..." for the
        # "gdc-client_2.3_Ubuntu_x64-....zip" release
        build = zip_name.split(f"_{label}_", 1)[1].removesuffix(".zip")
        name_version = name.split("_")[1]
        if (
            f"_{label}_" in name
            and name.endswith(build)
            and name_version.startswith(version)
        ):
            md5 = checksum
    return version, zip_name, md5


def install_gdc_client(version="2.3", install_dir="tools/gdc-client"):
    """Return gdc-client from ``{install_dir}/{version}/``, downloading it if missing.

    If the executable is already in the version folder (installed earlier, or
    put there by hand), it is used after checking that ``--version`` matches
    the folder. Otherwise the version folder is created and the release zip is
    downloaded from the GDC, md5-checked and unpacked.

    Parameters
    ----------
    version : str, optional
        A key of :data:`PINNED_RELEASES` (e.g. ``"2.3"``) or ``"latest"``.
    install_dir : str, optional
        Parent directory holding one folder per version.

    Returns
    -------
    dict
        ``path``, ``version`` (from ``--version``), ``source``
        (``"downloaded"`` or ``"existing"``), ``source_url``, ``zip_md5``,
        ``zip_md5_verified``, ``binary_md5``, ``binary_md5_matches_release``
        and ``recorded_on``. Also written to
        ``{install_dir}/{version}/gdc_client_info.json``.

    Raises
    ------
    ValueError
        If ``version`` is not pinned and not ``"latest"``, or the platform is
        not supported.
    RuntimeError
        If the downloaded zip's md5 differs from the published one, or the
        executable in the version folder reports a different version.
    """
    system = platform.system()
    if system not in GDC_PLATFORM:
        raise ValueError(f"gdc-client is not distributed for {system}")

    if version == "latest":
        version, zip_name, zip_md5 = latest_gdc_client_release(system)
        binary_md5 = PINNED_RELEASES.get(version, {}).get(system, (None, None, None))[2]
    elif version in PINNED_RELEASES:
        zip_name, zip_md5, binary_md5 = PINNED_RELEASES[version][system]
    else:
        raise ValueError(
            f"gdc-client version {version!r} is not pinned; use one of "
            f"{list(PINNED_RELEASES)} or 'latest'"
        )

    target_dir = os.path.join(install_dir, version)
    exe = "gdc-client.exe" if system == "Windows" else "gdc-client"
    path = os.path.join(target_dir, exe)
    info_path = os.path.join(target_dir, "gdc_client_info.json")
    url = f"{GDC_FILE_URL}/{zip_name}"

    if os.path.isfile(path):
        if os.path.isfile(info_path):
            with open(info_path) as f:
                return json.load(f)
        source, observed_zip_md5 = "existing", None
    else:
        os.makedirs(target_dir, exist_ok=True)
        zip_path = os.path.join(target_dir, zip_name)
        logger.info(
            f"gdc-client {version} not found in {target_dir}; downloading {url}"
        )
        urllib.request.urlretrieve(url, zip_path)

        observed_zip_md5 = file_md5(zip_path)
        if zip_md5 is not None and observed_zip_md5 != zip_md5:
            raise RuntimeError(
                f"md5 mismatch for {zip_name}: "
                f"expected {zip_md5}, got {observed_zip_md5}"
            )
        if zip_md5 is None:
            logger.warning(f"  no published md5 found for {zip_name}; not verified")

        # The GDC zip contains a second zip, which contains the executable
        with zipfile.ZipFile(zip_path) as outer:
            outer.extractall(target_dir)
            inner_names = [n for n in outer.namelist() if n.endswith(".zip")]
        for inner in inner_names:
            with zipfile.ZipFile(os.path.join(target_dir, inner)) as z:
                z.extractall(target_dir)
            os.remove(os.path.join(target_dir, inner))
        os.remove(zip_path)
        os.chmod(
            path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
        )
        source = "downloaded"

    found_version = gdc_client_version(path)
    if found_version != version:
        raise RuntimeError(
            f"{path} reports version {found_version}, expected {version}"
        )
    observed_binary_md5 = file_md5(path)
    info = dict(
        path=path,
        version=found_version,
        source=source,
        source_url=url,
        zip_md5=observed_zip_md5,
        zip_md5_verified=observed_zip_md5 is not None and zip_md5 is not None,
        binary_md5=observed_binary_md5,
        binary_md5_matches_release=None
        if binary_md5 is None
        else observed_binary_md5 == binary_md5,
        recorded_on=date.today().isoformat(),
    )
    with open(info_path, "w") as f:
        json.dump(info, f, indent=2)
    return info


def ensure_gdc_client(
    version="2.3", install_dir="tools/gdc-client", gdc_client_path=None
):
    """Return a working gdc-client.

    By default the client lives in ``{install_dir}/{version}/`` and is
    downloaded there if missing (:func:`install_gdc_client`). A different
    executable can be given with ``gdc_client_path``; if its ``--version``
    differs from ``version`` (and ``version`` is not ``"latest"``), a warning
    is printed and it is still used.

    Parameters
    ----------
    version : str, optional
        Pinned version (e.g. ``"2.3"``) or ``"latest"``.
    install_dir : str, optional
        Parent directory holding one folder per version.
    gdc_client_path : str or None, optional
        Optional path to a gdc-client outside ``install_dir``.

    Returns
    -------
    dict
        ``path``, ``version``, ``source`` (``"downloaded"``, ``"existing"`` or
        ``"user-provided"``), ``binary_md5`` and, for clients in
        ``install_dir``, the release URL and md5 checks.
    """
    if gdc_client_path is not None:
        if not (
            os.path.isfile(gdc_client_path) and os.access(gdc_client_path, os.X_OK)
        ):
            raise FileNotFoundError(
                f"gdc-client not found or not executable: {gdc_client_path}"
            )
        found = gdc_client_version(gdc_client_path)
        info = dict(
            path=gdc_client_path,
            version=found,
            source="user-provided",
            binary_md5=file_md5(gdc_client_path),
        )
        if version != "latest" and found != version:
            logger.warning(
                f"user-provided gdc-client is version {found}, "
                f"pinned version is {version}"
            )
    else:
        info = install_gdc_client(version=version, install_dir=install_dir)
    logger.info(
        f"Using gdc-client {info['version']} ({info['source']}): {info['path']}"
    )
    return info
