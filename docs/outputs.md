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

