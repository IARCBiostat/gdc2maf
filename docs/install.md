# Installation

## Requirements

Python 3.10 or newer, plus pandas, numpy and requests. No GDC token is needed: gdc2maf
handles open-access files only.

## pip

```bash
git clone https://github.com/IARCBiostat/gdc2maf.git
cd gdc2maf
pip install -e .
```

### For developers
Add the extras you need:

```bash
pip install -e '.[dev]'     # pytest
pip install -e '.[docs]'    # sphinx and this site's theme
```

`requirements.txt` lists the same runtime dependencies for environments that install them
separately, such as a Docker layer or a locked CI job:

```bash
pip install -r requirements.txt
pip install -e . --no-deps
```

## conda

```bash
conda env create -f environment.yml
conda activate gdc2maf
pip install -e .
```

## Platforms

Everything except the download is plain Python and platform-independent: the GDC queries,
the file selection, the merge and the QC run anywhere gdc2maf installs.

The download step shells out to the GDC Data Transfer Tool (`gdc-client`).
The GDC publishes `gdc-client` as **compiled bundles for Linux, macOS and Windows** — see
{data}`gdc2maf.client.PINNED_RELEASES`.
gdc2maf picks the build matching `platform.system()` and downloads it.