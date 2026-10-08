# Outputs and reruns


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
| `<name>_gdc_annotations.tsv` | Every GDC annotation on the cohort's files, aliquots and cases |
| `<name>_quality_flags.tsv` | Patients flagged by GDC curation or the PanCanAtlas `Do_not_use` flag, with the reason |
| `<name>_<spec>.maf` | The merged MAF, with the flagged patients left out unless `--keep-flagged` |
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

The name carries the cohort, so a download directory shared by several cohorts (see
below) collects one record per cohort side by side rather than overwriting anything —
which makes that directory the one place to read what every cohort in the study used.

Nothing is deleted or overwritten silently except the log: a rerun reuses the tables it
finds, and `--refresh` / `refresh=[...]` is how you ask for a step to be redone.

(one-download-dir)=
## One download directory for every cohort

`--download-dir` (`download_dir=`) is not per cohort, and should not be made per
cohort. gdc-client stores each file as **`DIR/<file_id>/<file_name>`** — named
after the GDC file UUID, with nothing about a cohort in the path — and the merge
looks its files up by UUID rather than by reading the directory. So one directory
serves any number of cohorts:

```bash
gdc2maf maf --project TCGA-LUAD TCGA-LUSC --name Lung_MALE  --download-dir data/mafs
gdc2maf maf --project TCGA-BRCA            --name Breast_F  --download-dir data/mafs
```

Each run downloads only the files its own cohort selected, and skips any that are
already there and pass the size and md5 check.

:::{warning}
Give each cohort its own download directory and you get **duplicate copies of the
same file** — every cohort that selects a file downloads it again, under its own
directory, and the GDC is asked for bytes you already have.

It costs nothing while cohorts are disjoint, which is why the mistake is easy to
miss: two sex-split or different-project cohorts select different files, so
nothing is duplicated and the layout looks fine. It bites as soon as one cohort
overlaps another — a sensitivity subset (one analyte only, or `exclude_flagged=False`
against the same patients), a combined cohort covering two existing ones, a
re-slicing of the same projects. Those select files you have already downloaded,
and a per-cohort directory fetches every one again.

A shared directory also re-verifies every file's md5 on each run instead of one
cohort's slice, so corruption in a file two cohorts share is caught either way.
:::

Moving to a shared directory needs no re-download. The `<file_id>/` folders are
self-contained — each holds its `.maf.gz` and gdc-client's own `logs/` — so they
can be moved up a level and the next run finds them by UUID and verifies them:

```bash
mkdir -p data/mafs
mv data/mafs_per_cohort/*/*/ data/mafs/
```

## What happens on a rerun

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

One rerun does need `--refresh merge`: turning the quality-flag exclusion on (or off)
after a cohort has already been merged. The merged MAF is cached, so it would otherwise
be reused as it was built, with the excluded patients still in it. The run says so:

```text
WARNING  Lung_MALE: out/Lung_MALE/Lung_MALE_wxs_ensemble.maf was merged before these
         patients were excluded and still contains 55 of them. Rerun with
         refresh=['merge'] (--refresh merge) to rebuild it, or set
         exclude_flagged=False (--keep-flagged).
```

