# Your first cohort

This builds one verified, merged MAF for male lung-cancer patients in TCGA, and
reads off what happened to every patient on the way. It takes a few minutes and
about 70 MB of disk.

## The short version

Male lung-cancer patients in TCGA are two projects, `TCGA-LUAD` and `TCGA-LUSC`
({ref}`finding-the-project-ids` shows how to pick them). The whole run is one call:

```python
from gdc2maf import Cohort, configure_logging

configure_logging()      # the package logs nothing until you ask

cohort = Cohort(
    name="Lung_MALE",
    projects=["TCGA-LUAD", "TCGA-LUSC"],
    sex_at_birth=["male"],
    out_dir="out/Lung_MALE",
    download_dir="downloads",
)

result = cohort.run()
print(result["attrition"].to_string(index=False))
```

That resolves the cohort, lists its MAF files, keeps one per patient, installs an
md5-checked `gdc-client`, checks the quality flags, downloads and verifies every
file, merges them with QC, and writes every report and the plain-text download
record. For this cohort it takes 1089 project cases down to an analysis set of
537, with a reason recorded for each of the 552 that did not make it.

`result` is a dict holding every table — `cases`, `maf_files`, `selected`,
`dropped`, `download_check`, `sample_qc`, `attrition`, `summaries` — so nothing
is hidden behind the single call. {doc}`../outputs` lists the files it wrote.

The equivalent on the command line:

```bash
gdc2maf maf --project TCGA-LUAD TCGA-LUSC --sex-at-birth male \
            --name Lung_MALE --out out --download-dir downloads
```

## Step by step

The rest of this page does exactly the same work one step at a time, which is how
to look at what each step produced and check it before going on. Every method
caches to `out_dir`, so running them in sequence costs no more than the single
call above.

### 1. Find the projects

Where `["TCGA-LUAD", "TCGA-LUSC"]` came from. `projects` wants GDC project IDs,
not a tissue name, so start from the site:

```bash
curl -G https://api.gdc.cancer.gov/projects \
  --data-urlencode 'fields=project_id,name' \
  --data-urlencode 'size=100' \
  --data-urlencode 'format=TSV' \
  --data-urlencode 'filters={"op":"in","content":{"field":"primary_site","value":["Bronchus and lung"]}}'
```

That returns 28 projects, because a *primary site* is not a disease: `TCGA-MESO`
(mesothelioma) and `TCGA-SKCM` (melanoma) are in the list too, each holding some
case recorded at a lung site. The lung-cancer projects are `TCGA-LUAD` and
`TCGA-LUSC`. {ref}`finding-the-project-ids` covers this in full.

### 2. Resolve the cohort

With the same `cohort` object built at the top of this page, ask for the cases
instead of calling `run()`:

```python
cases = cohort.cases()
print(f"{cases['included'].sum()} of {len(cases)} project cases are in the cohort")
```

Every case of both projects stays in the table. The ones that fail a criterion
are *marked*, not dropped:

```python
print(cases.loc[~cases["included"], ["submitter_id", "exclusion_step", "exclusion_reason"]].head())
```

```text
    submitter_id    exclusion_step          exclusion_reason
0   TCGA-05-4244    2. sex at birth         sex_at_birth is female
1   TCGA-05-4249    2. sex at birth         sex_at_birth is female
...
```

That is what makes the attrition report at the end possible.

### 3. List the candidate files

Nothing is downloaded yet:

```python
maf_files = cohort.maf_files()
print(len(maf_files), "files for", maf_files["case_id"].nunique(), "patients")
```

Some patients have several files — different tumour aliquots, different normals,
different sequencing centres. To see how they differ:

```python
duplicates = cohort.duplicates()
print(duplicates["pattern"].value_counts().head())
```

```text
same_tumour_sample_diff_aliquot + tumour_wga_and_native                            36
tumour_sample_differs + same_tumour_sample_diff_aliquot + tumour_wga_and_native     4
...
```

### 4. Keep one file per patient

```python
selected, dropped = cohort.selection()
print(dropped[["submitter_id", "drop_reason"]].head())
```

```text
    submitter_id    drop_reason
0   TCGA-38-4631    WGA tumour DNA (native DNA available)
1   TCGA-44-2656    WGA tumour DNA (native DNA available)
...
```

Every dropped file records *which* rule decided against it, so the choice is
auditable rather than arbitrary. The rule is documented under
{func}`gdc2maf.selection.select_one_file_per_case`.

### 5. Download and verify

The first call installs an md5-checked `gdc-client` if there is not one already:

```python
check = cohort.download(n_clients=8)
print(check["status"].value_counts())
```

```text
ok    596
```

Every file is checked by size and md5 against the GDC's own metadata, before and
after the transfer. A rerun re-checks and downloads only what is missing or
corrupt.

Files land in `downloads/<file_id>/<file_name>` — named after the GDC file UUID,
with nothing about the cohort in the path. That is why `download_dir="downloads"`
above is **not** `"downloads/Lung_MALE"`: point every cohort you ever build at the
same directory and each file is downloaded once, however many cohorts select it.
The merge finds its own files by UUID, so a directory holding other cohorts'
files costs nothing.

```python
Cohort(name="Lung_MALE",  ..., download_dir="downloads")
Cohort(name="Breast_F",   ..., download_dir="downloads")   # the same directory
```

A directory per cohort instead gives you duplicate copies of the same file as
soon as two cohorts overlap — see {ref}`one-download-dir`.

### 6. Check the quality flags

Two lists say a patient's material may be unusable: the GDC's own curation
annotations, and the PanCanAtlas `Do_not_use` flag. Both are consulted without
any setup — the annotations come from the API, and the PanCanAtlas table is
downloaded once into `pancan_dir`, md5-checked:

```python
print(cohort.quality_flags())
```

```text
    submitter_id       source      reason
0  TCGA-05-4382  PanCanAtlas  Do_not_use
1  TCGA-05-4418  PanCanAtlas  Do_not_use
...
```

56 of the 596 selected patients are flagged here, all by the PanCanAtlas. They
are excluded by default: the next step leaves their files out of the merge, so
they are absent from the MAF and from everything built on it. The table is saved
as `Lung_MALE_quality_flags.tsv` either way.

`--keep-flagged` (or `exclude_flagged=False`) keeps them and only reports the
flags. {doc}`../quality` has the evidence behind the 56 and the case for each
choice — worth reading once before you rely on either.

### 7. Merge, with QC

```python
maf, sample_qc, qc_summary = cohort.merged_maf()
print(len(maf), "variants from", sample_qc["has_variants"].sum(), "patients")
```

The merge is strict about anything that would corrupt the result and raises
rather than guessing — a failed download check, a header that does not match, a
row on the wrong reference build, a barcode that is not the selected aliquot, a
duplicated variant.

It is deliberately *not* strict about biology. Samples with an implausible
mutation count, excess low-VAF calls or weak caller support are flagged with the
evidence and left in the MAF:

```python
print(cohort.samples_for_review())
```

Whether to drop them is yours to decide, not the loader's.

### 8. Read the attrition

```python
print(cohort.attrition().drop(columns="cohort").to_string(index=False))
```

```text
                          step                                          reason  n_removed  n_remaining
             1. GDC project(s)                                       all cases          0         1089
               2. sex at birth                          sex_at_birth is female        411          678
               2. sex at birth                sex_at_birth missing in GDC         63          615
      3. open WXS ensemble MAF                        no open WXS ensemble MAF         19          596
           4. MAF file listing                                               -          0          596
5. download (size + md5 check)                                               -          0          596
              6. quality flags                          flagged by PanCanAtlas         56          540
                 7. merged MAF         MAF has no variants (after GDC masking)          3          537
```

Read it top to bottom: the cohort starts as every case of both projects and
narrows one criterion at a time, each row naming what was lost and why. Steps
that removed nothing are still shown, so the arithmetic is checkable rather than
implied.

Each patient is counted once, at the step that lost them: four patients have no
variants after GDC masking, but one of them was already gone at step 6, so the
merge row removes three. With `--keep-flagged` step 6 disappears and the merge
row removes all four, ending at 592.

Counts are from GDC Data Release 46.0; a later release will differ. Had you
passed `sex_fallback_path` ({ref}`sex-criterion`), cases with no sex in the GDC
but a sex in that table would be filled in before this step — for this cohort
there are none, so the numbers would be the same.

## Next

- {doc}`signature-pipeline` builds a complete analysis on top of this.
- {doc}`../outputs` lists every file the run wrote.
- {doc}`../cohorts` covers cohort criteria and non-WXS file kinds.
