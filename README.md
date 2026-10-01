# gdc2maf

From a GDC cohort specification to one verified, merged MAF.

Say which GDC project(s) and which cases you want. `gdc2maf` resolves the cohort against
the GDC API, lists the matching open-access MAF files, picks one file per patient,
installs an md5-verified `gdc-client`, downloads and checks every file, and merges them
into a single MAF with per-sample QC — accounting for every patient lost along the way,
with the reason.

## Quick start

Clone the project and install:
```bash
pip install -e .
```

To resolve and download your cohort from the GDC's open-access MAF files (token based download will be supported in the future) use eighther of the two methods below:

**Command line** — one command, start to finish:

```bash
gdc2maf maf --project TCGA-LUAD --sex-at-birth male --name Lung_MALE
```

Writes `out/Lung_MALE/` (cases, file listing, selection, merged MAF, QC and an attrition
table) and downloads into `downloads/`. See [Command line](#command-line) for the
per-step subcommands and the full option list.

**Python** — the same run, as an object:

```python
from gdc2maf import Cohort, configure_logging

configure_logging()
cohort = Cohort(
    name="Lung_MALE",
    projects=["TCGA-LUAD"],
    sex_at_birth=["male"],
    out_dir="out/Lung_MALE",
    download_dir="downloads",
)
result = cohort.run()
print(result["attrition"])
```

Every step is also a method that returns DataFrames, so you can stop anywhere and inspect
what happened. See [Python](#python) for the step-by-step form.

# Detailed examples

## Command line

The quick-start run with every path spelled out — two projects, explicit output and
download directories, eight parallel download processes:

```bash
gdc2maf maf --project TCGA-LUAD --project TCGA-LUSC \
            --sex-at-birth male --name Lung_MALE \
            --out out --download-dir downloads --jobs 8
```

That writes everything to `out/Lung_MALE/`, downloads into `downloads/<file_id>/`, and
prints the attrition table and a summary of each step.

**The `gdc-client` is handled for you.** Downloading GDC files needs the GDC Data Transfer
Tool, so before the first download gdc2maf looks for `gdc-client` in
`tools/gdc-client/<version>/` and, if it is not there, downloads that release from the GDC
and checks its md5 against a pinned table before running it. Which binary ran — path,
version, and its own md5 — is recorded in `<name>_gdc_client.json` beside the cohort's
other outputs. Version `2.3` is pinned by default; `--gdc-client-version latest` takes the
newest release, `--gdc-client PATH` uses a client you already have, and `--install-dir`
moves where versions are kept. `gdc2maf client` does only this step, which is useful for
warming a shared installation before a batch of cluster jobs.

A rerun reuses what is already on disk, so it is cheap to run again after an interruption:

```bash
gdc2maf maf ... --refresh files      # re-query the file listing, keep the rest
gdc2maf maf ... --no-download        # verify the files on disk, download nothing
```

Each step is also its own subcommand, for splitting a run across cluster jobs:

| Command | Does |
|---|---|
| `gdc2maf cases` | Resolve the cohort into a table of GDC cases |
| `gdc2maf clinical` | Fetch clinical data for those cases (auxiliary) |
| `gdc2maf files` | List the MAF files; describe patients with more than one |
| `gdc2maf select` | Choose one file per patient, with the reason each other was dropped |
| `gdc2maf download` | Download the selected files; verify size and md5 |
| `gdc2maf merge` | Merge the verified MAFs into one MAF, with QC |
| `gdc2maf reports` | Rebuild the attrition and summary tables from cached outputs |
| `gdc2maf client` | Locate or install `gdc-client` and report which binary is used |
| `gdc2maf maf` | All of the above, in order |

Every command has built-in help: `gdc2maf -h` (or `--help`) lists the commands, and
`gdc2maf <command> -h` lists that command's options with their defaults. Useful ones: `--sex-fallback` (fill
`sex_at_birth` from an external clinical table where the GDC has none),
`--remove-gdc-filter-flag` (drop `GDC_FILTER`-flagged variants; the default keeps them and
reports the counts), `--gdc-client-version` (pin a version, or `latest`), `--no-log`,
`--quiet`.

### Finding the project IDs

`--project` takes GDC project IDs such as `TCGA-LUAD` — not a tissue or a disease name.
Two ways to look one up, both repeated in `gdc2maf -h`:

1. Browse <https://portal.gdc.cancer.gov/projects> and filter by primary site or program;
   the **Project** column holds the IDs.
2. Ask the API. Every project holding a lung case, as a table:

   ```bash
   curl -G https://api.gdc.cancer.gov/projects \
     --data-urlencode 'fields=project_id,name' \
     --data-urlencode 'size=100' \
     --data-urlencode 'format=TSV' \
     --data-urlencode 'filters={"op":"in","content":{"field":"primary_site","value":["Bronchus and lung"]}}'
   ```

   Drop the last line to list every project. Site names follow ICD-O-3, so lung is
   `Bronchus and lung`, not `Lung`.

**A primary site is not a disease.** That lung query returns 28 projects — `TCGA-MESO`
(mesothelioma) and `TCGA-SKCM` (melanoma) among them — because each holds at least one case
recorded at a lung site. This is why gdc2maf asks for projects rather than a site: pick the
projects for the disease you mean and pass each with its own `--project`. TCGA lung cancer
is `TCGA-LUAD` (adenocarcinoma) plus `TCGA-LUSC` (squamous cell):

```bash
gdc2maf maf --project TCGA-LUAD --project TCGA-LUSC --name Lung
```

Every case of the named projects is fetched and kept in the cases table, so the attrition
report starts from the whole project and shows what each criterion removed.

### Choosing which cases

`--sex-at-birth` is optional: **omit it and every sex is kept**, including cases whose sex
the GDC does not record. Give it once (`--sex-at-birth male`) or repeat it to restrict the
cohort; cases of any other sex, and cases with no recorded sex, are then excluded and
counted in the attrition table with the reason. The Python equivalent is
`sex_at_birth=None`, which is also the default.

## Python

Step by step instead of `cohort.run()`. Each method returns DataFrames and caches to
`out_dir`; calling a late step runs the earlier ones it needs, so you can start anywhere.

```python
from gdc2maf import Cohort, configure_logging

configure_logging()  # the package logs nothing until you ask

cohort = Cohort(
    name="Lung_MALE",
    projects=["TCGA-LUAD", "TCGA-LUSC"],
    sex_at_birth=["male"],
    out_dir="out/Lung_MALE",
    download_dir="downloads",
)

cases = cohort.cases()  # GDC query, or the cached TSV
maf_files = cohort.maf_files()  # one row per candidate file, nothing downloaded
selected, dropped = cohort.selection()  # one file per patient + why the others went
check = cohort.download(n_clients=8)  # gdc-client, then size + md5 on every file
maf, sample_qc, qc = cohort.merged_maf()

print(cohort.attrition().to_string(index=False))
```

`cohort.run()` does all of it and returns every table in a dict. To work from a merged MAF
that is already on disk without loading all of it, use
`cohort.read_merged_maf(usecols=[...])`.

The step functions are importable on their own (`fetch_cases`, `list_maf_files`,
`select_one_file_per_case`, `download_files`, `merge_mafs`, `ensure_gdc_client`, …) if you
would rather drive them yourself than hold a `Cohort`.

## Which files are "the cohort's MAFs"

By default, the open-access GRCh38 `Masked Somatic Mutation` files from the GDC's Aliquot
Ensemble Somatic Variant Merging and Masking workflow, exome (WXS) — one MAF per
tumour-normal aliquot pair. That is a value, not hard-coded logic:

```python
from gdc2maf import Cohort, FileSpec

wgs_ensemble = FileSpec(
    name="wgs_ensemble_maf",
    title="open WGS ensemble MAF",
    data_category=("Simple Nucleotide Variation",),
    data_type=("Masked Somatic Mutation",),
    experimental_strategy=("WGS",),
    access=("open",),
)
cohort = Cohort(..., spec=wgs_ensemble)
```

Only `WXS_ENSEMBLE_MAF` is exercised by this package's own runs; other specs are passed
straight through to the GDC `/files` endpoint.

## Outputs

Everything lands in `out_dir`, prefixed with the cohort name:

| File | Contents |
|---|---|
| `<name>_cases.tsv` | Every case of the project(s), with `exclusion_step`, `exclusion_reason`, `included` |
| `<name>_query.json` | The exact query: projects, criteria, file spec, GDC data release, date |
| `<name>_clinical.tsv` | Demographic / diagnosis / follow-up / exposure data per case |
| `<name>_files_metadata.tsv` | Every candidate MAF file, with md5 and tumour/normal aliquot barcodes |
| `<name>_duplicate_patients.tsv` | Patients with >1 file, and how the files differ |
| `<name>_files_selected.tsv`, `<name>_files_dropped.tsv` | The chosen file per patient, and why each other was dropped |
| `<name>_manifest.txt` | gdc-client manifest of the selected files |
| `<name>_download_check.tsv` | Size and md5 status of every selected file |
| `<name>_gdc_client.json` | Which gdc-client binary ran, its version and md5 |
| `<name>_<spec>.maf` | The merged MAF |
| `<name>_sample_qc.tsv` | Per-sample metrics and artefact flags |
| `<name>_samples_for_review.tsv` | Samples flagged for manual review, with the evidence |
| `<name>_case_attrition.tsv` | Cases lost per step, by reason |
| `<name>_*_summary.tsv` | One-row summaries, stackable across cohorts with `summary_table` |
| `<name>_run.log` | Everything logged during the run (overwritten each run) |

And one file goes next to the **downloaded** MAFs, in `--download-dir`, because that
directory is otherwise just numbered folders of `.maf.gz` with nothing to say what they
are:

| File | Contents |
|---|---|
| `<name>_download_record.txt` | Plain-text record: the cohort, the GDC release, the gdc-client used, how many cases were excluded at each step and why, how many files were downloaded and verified, and how many patients ended up with variants |

Pass `--record-dir` (or `record_dir=` to `cohort.run()`) to put it somewhere else, such as
beside the merged MAF in `out_dir`.

Nothing is deleted or overwritten silently except the log: a rerun reuses the tables it
finds, and `--refresh` / `refresh=[...]` is how you ask for a step to be redone.

### What happens on a rerun

Running the same command again does not redo the work:

- **Each step reuses its cached table** if the file is already in `out_dir`, so no GDC
  query is repeated. `--refresh files` (or `refresh=["files"]`) redoes one step;
  `--refresh all` redoes all of them.
- **Files already on disk are not downloaded again.** Every selected file is checked by
  size and md5 first, and only the missing or corrupt ones are fetched.
- **The download record is left untouched** when nothing changed. It holds no timestamp,
  so its text is identical whenever the result is, and the writer compares before writing:

  ```python
  cohort.run()["record"]    # {'status': 'created',   'superseded': None}
  cohort.run()["record"]    # {'status': 'unchanged', 'superseded': None}
  ```

  When the content does change — a new GDC release, a widened cohort, a file that fails
  verification — the old record is first renamed to
  `<name>_download_record_superseded_<date>.txt`, so nothing is lost.
- **Only `<name>_run.log` is overwritten**, since it describes the current run.

## Notes

- **Nothing is dropped quietly.** Cases that fail a criterion stay in the cases table,
  marked with the step and reason, so the attrition table runs from all project cases down
  to the merged MAF. Dropped duplicate files record which selection key decided against
  them.
- **Strict about corruption, permissive about biology.** The merge raises on a failed
  download check, a mismatched header, a non-GRCh38 row, a barcode that is not the
  selected aliquot, or a duplicated variant. Samples with implausible counts, excess
  low-VAF calls or weak caller support are flagged with the evidence and left in — whether
  to drop them is the analyst's call.
- **Verified binaries.** `gdc-client` releases are md5-checked against a pinned table
  before use, and the binary's own md5 is recorded with the run.
- **Reproducible queries.** `<name>_query.json` records the GDC data release and the exact
  cohort criteria, so a table can be traced back to the release that produced it.
- **Logging, not printing.** The package logs to the `gdc2maf` logger and configures
  nothing on import. `configure_logging()` attaches handlers; `tee_stdout()` is there for
  when a downstream tool's own stdout has to reach the log file too.
