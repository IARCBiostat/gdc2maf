# Why it works this way

Every step of the cohort-to-MAF path had alternatives. This records what was
chosen, what was not, and why — so that a result can be defended, and so that
changing a choice is a decision rather than an accident.

Two principles run through all of it. **Nothing is dropped quietly**: a case
that fails a criterion is marked, not deleted, so the attrition table runs from
every project case down to the analysis set with a reason for each loss. And
**no choice looks at the mutations**: a rule like "keep the file with the most
variants" would inflate tumour mutational burden and bias the signatures it is
meant to measure.

## Defining the cohort

**Projects, not primary site.** A cohort is GDC project IDs. Filtering on
`primary_site` instead sweeps in cases from projects that merely contain a case
biopsied at that site — a lung-site query returns mesothelioma and melanoma
projects alongside the lung-cancer ones. See {ref}`finding-the-project-ids`.

**An API query, not a portal export.** The query is built from parameters and
recorded in `<name>_query.json` with the GDC data release, so a table can be
traced to what produced it. A portal export is a manual step whose filters live
only in a browser session.

**`demographic.sex_at_birth`, filtered in pandas.** `demographic.gender` is no
longer in the GDC data dictionary and returns nothing. The filter is applied
after the query rather than inside it, so excluded cases stay in the table and
can be counted — the attrition principle above.

**A missing sex can be filled from an external table.** Three options were
considered:

1. *Chosen:* fill from a clinical table such as the PanCanAtlas one
   ({func}`gdc2maf.cases.fill_sex_from_table`), with the GDC value always
   winning, conflicts logged, and the source recorded in
   `sex_at_birth_source`. Cases still missing a sex are excluded with the
   reason recorded.
2. Exclude every case with no GDC value, with no fallback. Available:
   leave `sex_fallback_path` unset.
3. Infer genetic sex from sequencing data, e.g. chrY coverage. Not done — it
   needs controlled-access BAMs, and genetic sex is not the recorded clinical
   variable a cohort is defined on.

## Which files are the cohort's MAFs

**Six conditions, all on the same file** — see
{data}`gdc2maf.spec.WXS_ENSEMBLE_MAF`. These are the open-access GRCh38
calls from the [Aliquot Ensemble Somatic Variant Merging and
Masking](https://docs.gdc.cancer.gov/Encyclopedia/pages/Aliquot_Ensemble_Somatic_Variant_Merging_and_Masking/)
workflow, which merges the per-caller MAFs (MuSE, MuTect2, VarScan2, Pindel)
and removes potentially identifiable germline variants. The GDC notes that
those open-access criteria "over-filter some of the true positive somatic
variants" to keep germline information private.

Not taken:

- **Per-caller VCFs** (`Raw Simple Somatic Mutation`, `Annotated Somatic
  Mutation`): single callers, and controlled access.
- **The portal cohort MAF**: the same calls merged for a portal cohort, but you
  cannot control which aliquot is used per patient, and it is a manual
  download.
- **MC3** (`mc3.v0.2.8.PUBLIC.maf`): the earlier GRCh37 call set.
- **WGS**: a different dataset. An exome analysis must also run SigProfiler
  with `exome=True`.

Any of these is a {class}`~gdc2maf.spec.FileSpec` away; see
{doc}`cohorts`.

**`/files`, not `/cases`, for "does this case have a MAF".** On `/cases`,
nested `files.*` conditions can be satisfied by *different* files of the same
case — one WXS file and a separate masked MAF would pass a combined filter.
`/files` applies every condition to one file.
{func}`gdc2maf.api.cases_with_files` therefore queries `/files` and reads the
case IDs back off the hits.

**Tumour and normal from the GDC's own annotation.** Each file links to exactly
two aliquots; each aliquot maps to a sample, and the sample's `tissue_type`
gives the role. Parsing the barcode's sample-type code would work for TCGA, but
the GDC annotation is direct and also records the sample type label — `FFPE
Scrolls` carries code `01` like any other primary tumour. The code raises
unless a file has exactly one case, one tumour and one normal.

## One file per patient

A patient sequenced more than once has more than one ensemble MAF, and
signature analysis needs exactly one sample per patient. The rule has to be
fixed in advance, for the reason given at the top: anything that reads the
mutations biases the result.

{func}`gdc2maf.selection.select_one_file_per_case` sorts a patient's files by
these keys and keeps the first:

| # | Prefer | Over | Why |
|---|---|---|---|
| 1 | Primary tumour (`01`, `03`, `09`) | Metastatic, recurrent, new primary | Matches the primary-tumour question and TCGA convention |
| 2 | Non-FFPE tumour | FFPE tumour | Formalin fixation causes C>T artefacts that look like SBS signatures |
| 3 | Native DNA (`D`) | WGA DNA (`W`, `X`) | Amplification causes artefacts and uneven coverage |
| 4 | Lower vial (`A`) | Higher vial (`B`, `C`) | TCGA convention; fixed and arbitrary |
| 5 | Later plate | Earlier plate | Later plates are usually re-sequencing with newer chemistry; plate codes compare as text (`A27T` > `A271` > `1753`) |
| 6 | Lower portion | Higher portion | Fixed and arbitrary |
| 7 | Blood-derived normal (`10`) | Solid-tissue normal (`11`) | Adjacent normal tissue can contain tumour cells, which masks true somatic variants |
| 8 | Tumour, then normal aliquot barcode, alphabetically | — | Final tie-break, so the result never depends on input row order |

The properties that matter:

- **Deterministic.** The keys give a total order, and the code raises if two of
  a patient's files share both aliquots. Shuffling the input does not change
  the selection — there is a test for exactly that.
- **Self-contained.** Only GDC metadata is used, so anyone can reproduce it.
- **Consistent with the GDC's own approach.** For its aggregated open-access
  MAFs the GDC also keeps "only one tumor-normal pair … based on the plate
  number, sample type, analyte type and other features extracted from tumor
  TCGA aliquot barcode" ([MAF
  format](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/)). The
  exact ordering is not published, so this rule is written out in full.
- **Auditable.** Every dropped file records which key decided against it, in
  `<name>_files_dropped.tsv`.

The keys read TCGA aliquot barcodes ({mod}`gdc2maf.tcga`). For a non-TCGA
project, sort `maf_files` yourself and keep one row per `case_id`.

## Downloading

**gdc-client, driven by a manifest.** It retries and can resume. Not taken:
`POST /data` with batches of IDs, which returns a tar.gz to unpack and match
back to IDs; and the portal's download button, which is manual. See
{doc}`install` on which platforms have a published build.

**Several gdc-client processes, not one.** One process fetches a manifest's
files one after another — a couple of seconds each for MAFs of a few tens of
KB, which runs to hours for a few thousand files. The tool's `-n` option
("number of client connections") did not help noticeably for many small files.
Splitting the manifest across 8 processes was roughly seven times faster on a
test batch. Lower `n_clients` if the GDC rate-limits.

**`--no-annotations --no-related-files`.** Otherwise gdc-client writes an
`annotations.txt` beside every MAF. Annotations are better queried from the API
when wanted — see {doc}`quality`.

**A size and md5 check of our own, every run.** gdc-client checks md5 during
transfer, but a separate check proves that the files on disk *now* are complete
and are the selected ones. Size is compared first and md5 computed only if it
matches. "Verified = selection" compares the *sets* of file IDs, not just the
counts. A non-zero gdc-client exit code is not fatal: the check decides, and
failures are counted in the attrition.

## Merging, and what counts as bad

The merge is strict about anything that would silently corrupt the result and
raises: a file that failed its download check, a header that differs between
files, a row on the wrong reference build, a barcode that is not the selected
aliquot, a duplicated variant within a sample.

It is deliberately *not* strict about biology. A high mutation count cannot by
itself distinguish an error from a real hypermutator — POLE, mismatch-repair
and APOBEC activity, and UV or tobacco exposure, all produce genuinely high
counts. So {func}`gdc2maf.maf.flag_samples` flags samples **with the evidence**
and leaves them in:

| Flag | Rule |
|---|---|
| `high_count` | log10(`n_variants`) > Q3 + 3 × IQR within the cohort — the upper outer fence, on the log scale |
| `low_vaf_excess` | ≥ 25% of variants with VAF < 0.1 |
| `weak_caller_support` | ≥ 35% of variants called by ≤ 2 callers |
| `review` | `high_count`, **or** top 5% of counts **and** (`low_vaf_excess` or `weak_caller_support`) |

alongside the substitution spectrum, which points at the likely process: C>A at
low VAF suggests an oxidative artefact, C>G with C>T suggests APOBEC, very high
C>T at high VAF suggests UV. Not taken: a fixed count cut-off, which would
remove real hypermutators and cannot hold across cohorts whose burden differs
several-fold; automatic removal of every `high_count` sample; and MC3-style
per-variant filters, which the GDC masked MAF does not carry.

Whole-genome amplification is reported in the `wga` column but does not trigger
review on its own, since it is a property of the library rather than evidence
of a bad call.

**`GDC_FILTER`-flagged variants are kept by default** and their per-flag counts
reported. The GDC's MAF format documentation lists `NonExonic` and `gdc_pon`
among variants removed when protected MAFs become open-access, yet release 46.0
ensemble MAFs still contain them; `common_in_gnomAD` marks possible germline
calls. Which of those to remove is a study's decision, so
`remove_gdc_filter_flags` takes the list and the default is to remove none.

## Quality flags

**Both lists are consulted, and acted on.** The GDC's curation annotations and
the PanCanAtlas `Do_not_use` flag are fetched without being asked — the
PanCanAtlas table is downloaded once, md5-checked — and the patients they flag
are left out of the merge by default. The alternative, reporting and keeping
them, was the first implementation and is still one switch away
(`--keep-flagged`); the default changed because a signature analysis is the
wrong place to discover that a sample was redacted or flagged do-not-use.

**Excluded at the merge, not in the report.** The filter is in
{meth}`~gdc2maf.cohort.Cohort.merged_maf`, so the patients are absent from the
MAF itself and from everything built on it. Doing it in the attrition report
instead would have been easier and wrong: the table would have said a patient
was dropped while their variants were still in the MAF. For the same reason the
exclusion is step 6, before the merge, rather than a step appended after it.

**Not in the selection, because it cannot be.**
{meth}`~gdc2maf.cohort.Cohort.annotations` queries the annotations of the
*selected* files, so filtering the selection on them would be circular. The
merge is the first point where every flag is known.

{doc}`quality` has the evidence behind the flags, what the default costs on a
real cohort, and why the mere presence of an annotation is not grounds for
exclusion.
