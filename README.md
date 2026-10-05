# gdc2maf

<!-- docs-include: intro-start -->
From a GDC cohort specification to one verified, merged MAF.

Say which GDC project(s) and which cases you want. `gdc2maf` resolves the cohort against
the GDC API, lists the matching open-access MAF files, picks one file per patient,
installs an md5-verified `gdc-client`, downloads and checks every file, and merges them
into a single MAF with per-sample QC — leaving out the patients the GDC's own curation or
the TCGA PanCanAtlas `Do_not_use` flag calls unusable, and accounting for every patient
lost along the way, with the reason. Downloading needs no GDC token: only open-access
files are handled.
<!-- docs-include: intro-end -->

**Full documentation — the option reference, the API reference and worked tutorials
including a complete mutational-signature pipeline — is at
<https://iarcbiostat.github.io/gdc2maf/>. Check the website for more information and
examples beyond the quick start below.**

<!-- docs-include: quickstart-start -->
## Quick start

Clone the project and install:

```bash
pip install -e .
```

To resolve and download your cohort from the GDC's open-access MAF files (token based
download will be supported in the future) use either of the two methods below:


**Python** — the same run, as an object:

```python
from gdc2maf import Cohort, configure_logging

configure_logging()
cohort = Cohort(
    name="Lung_MALE",
    projects=["TCGA-LUAD"],   # see Cohorts and file specs, linked above
    sex_at_birth=["male"],
    out_dir="out/Lung_MALE",
    download_dir="downloads",
)
result = cohort.run()
print(result["attrition"])
```

Every step is also a method that returns DataFrames, so you can stop anywhere and inspect
what happened.


**Command line** — one command, start to finish:

```bash
gdc2maf maf --project TCGA-LUAD --sex-at-birth male --name Lung_MALE
```

`--project` takes GDC project IDs, not a tissue or disease name, and a disease usually
spans more than one project — TCGA lung cancer is `TCGA-LUAD` *and* `TCGA-LUSC`.
[Cohorts and file specs](https://iarcbiostat.github.io/gdc2maf/cohorts.html#finding-the-project-ids)
shows how to look the right ones up on the GDC portal or from its API.

Writes `out/Lung_MALE/` (cases, file listing, selection, quality flags, merged MAF, QC
and an attrition table) and downloads into `downloads/`.
<!-- docs-include: quickstart-end -->

## Where to read more

Everything beyond the quick start lives in the documentation site, so it is written and
maintained in one place only:

| Topic | Page |
|---|---|
| Installing, and which platforms gdc-client supports | [Installation](https://iarcbiostat.github.io/gdc2maf/install.html) |
| Every command and option, and how to find project IDs | [Command line](https://iarcbiostat.github.io/gdc2maf/cli.html) |
| Driving the steps from Python | [Python API](https://iarcbiostat.github.io/gdc2maf/api.html) |
| Defining a cohort; which GDC files count as its MAFs | [Cohorts and file specs](https://iarcbiostat.github.io/gdc2maf/cohorts.html) |
| What each output file holds, and what a rerun does | [Outputs and reruns](https://iarcbiostat.github.io/gdc2maf/outputs.html) |
| Worked examples, including a full signature-extraction pipeline | [Tutorials](https://iarcbiostat.github.io/gdc2maf/tutorials/index.html) |
| Every class and function | [API reference](https://iarcbiostat.github.io/gdc2maf/reference/index.html) |

## Development

```bash
pip install -e '.[dev]'
pytest                      # unit tests; no network, no GDC access needed

pip install -e '.[docs]'
sphinx-build -b html docs docs/_build/html
```

## License

Copyright (C) 2026 Ali Farnudi.

gdc2maf is free software under the GNU General Public License, version 3 or later, and
comes with NO WARRANTY. See [LICENSE](LICENSE), or
<https://www.gnu.org/licenses/gpl-3.0.html>. The runtime dependencies (pandas and numpy
under BSD-3-Clause, requests under Apache-2.0) are all GPL-compatible.
