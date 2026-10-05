# Quality flags and exclusions

Two independent lists say a patient's material may be unusable: the GDC's own
curation annotations, and the TCGA PanCanAtlas `Do_not_use` flag. gdc2maf
fetches both, reports both, and by default excludes the patients either one
flags — their files are left out of the merged MAF, so nothing downstream sees
them.

Nothing has to be set up for that. The GDC annotations come from the API with
the rest of the cohort, and the PanCanAtlas table is downloaded on first use,
md5-checked, into `pancan_dir`:

```python
from gdc2maf import Cohort

cohort = Cohort(
    name="Lung_MALE",
    projects=["TCGA-LUAD", "TCGA-LUSC"],
    out_dir="out/Lung_MALE",
    download_dir="downloads",
)

result = cohort.run()             # flagged patients already excluded
flags = cohort.quality_flags()    # the record of who, and why
```

```bash
gdc2maf maf --project TCGA-LUAD TCGA-LUSC --name Lung_MALE
```

`flags` is one row per flag — `submitter_id`, `source`, `reason` — saved as
`<name>_quality_flags.tsv` whether or not the patients are excluded. The
exclusion itself shows up as step 6 of the attrition table, before the merge,
because that is where the files were left out:

```text
5. download (size + md5 check)                                       -     0    596
              6. quality flags                 flagged by PanCanAtlas      56    540
                 7. merged MAF    MAF has no variants (after GDC masking)   3    537
```

## The knobs

| | Python | Command line |
|---|---|---|
| Keep the flagged patients | `exclude_flagged=False` | `--keep-flagged` |
| Use a table already on disk | `quality_annotations_path=...` | `--quality-annotations TSV` |
| Where it is downloaded | `pancan_dir="data/pancan"` | `--pancan-dir DIR` |
| Do not use the PanCan flag at all | `fetch_quality_annotations=False` | `--no-quality-annotations` |

The download happens once and is reused from `pancan_dir` on every later run,
by every cohort pointed at the same directory. It is skipped altogether for a
cohort with no TCGA project — the PanCanAtlas covers TCGA only — and the run
log says so.

:::{warning}
The merge is cached. Turning the exclusion on or off after a cohort has already
been merged needs `--refresh merge` (`refresh=["merge"]`), or the MAF on disk is
reused exactly as it was built. A run that notices this says so:

```text
WARNING  Lung_MALE: ...Lung_MALE_wxs_ensemble.maf was merged before these
         patients were excluded and still contains 55 of them. Rerun with
         refresh=['merge'] (--refresh merge) to rebuild it, or set
         exclude_flagged=False (--keep-flagged).
```
:::

## Why the flagged patients go, and when to keep them

The default is the conservative reading — the same one the TCGA analysis
working groups applied — because a signature analysis is the wrong place to
discover that a sample was redacted or flagged do-not-use.

But keeping them is a defensible choice, because neither list is a verdict on
exome data.

The GDC's categories are mostly clinical or curatorial facts — `prior
malignancy`, `neoadjuvant therapy`, `item is noncanonical` — which a study may
well want to keep. The PanCanAtlas flag is per patient × aliquot × **platform**,
and most of those platforms are not exome sequencing; its reasons include prior
treatment, another platform's QC, duplicates and FFPE.

### What it costs, measured

On the male lung cohort of {doc}`tutorials/first-cohort`, 56 of 596 selected
patients are flagged, all by the PanCanAtlas — no selected file or aliquot
carried an exclusion-grade GDC annotation. Not one of the 56 is flagged on a
`WXS` row: `do_not_use_patients(annotations, platform="WXS")` is empty for the
whole table. The platforms behind them are

```text
Genome_Wide_SNP_6         32     HumanMethylation27          8
HumanMethylation450       13     IlluminaGA_miRNASeq         4
IlluminaHiSeq_RNASeqV2    13     IlluminaHiSeq_DNASeqC       2
IlluminaHiSeq_miRNASeq    12     MDA_RPPA_Core               9
```

so the evidence is mixed. Some of it is patient-level and applies however the
DNA was sequenced — "History of unacceptable prior treatment related to a
prior/other malignancy", an AWG pathology exclusion ("NSCLC, NOS — to be
excluded"), FFPE material. Some of it is about one aliquot on one platform: 32
of the 56 are there because a SNP6 copy-number aliquot failed Broad's pipeline
QC, which says nothing about the exome.

That is the trade the default makes: it accepts losing patients on non-exome
evidence rather than keeping one whose material the working groups rejected. If
your study wants the other side of it, `--keep-flagged` keeps them all and the
flags table still records every one.

So what a flag means for your analysis is a judgement about your analysis.
`--keep-flagged` makes it, and the flags table is written either way, so a
study that wants those patients back still has the full record of who was
flagged and on whose authority. What you should not do is leave the two out of
step: report an attrition table that drops them while analysing a MAF that
keeps them.

(annotation-presence)=
## An annotation is not an exclusion

Every open WXS ensemble MAF in the TCGA projects carries exactly one GDC
annotation, reading:

```text
Notification / General: Variants from SomaticSniper are not included.
```

That is a note about how the file was built. Treating the presence of an
annotation as grounds for exclusion would drop **every patient in the cohort** —
measured on the three documented cohorts, 954 of 958, 453 of 478 and 596 of 596
selected files carry one.

{func}`gdc2maf.annotations.exclusion_flags` is therefore the only thing that
decides, and it reads two narrow sets:

- {data}`gdc2maf.annotations.EXCLUSION_CLASSIFICATIONS` — `redaction` alone.
  It is the strongest signal the API carries; `notification`,
  `centernotification` and `observation` are commentary.
- {data}`gdc2maf.annotations.EXCLUSION_CATEGORIES` — the categories that
  question the identity or integrity of the material: `item flagged dnu`,
  `biospecimen identity unknown`, `genotype mismatch`, `sample compromised`,
  `duplicate item`.

Widen either per call if your study takes a different view:

```python
from gdc2maf import exclusion_flags

exclusion_flags(cohort.annotations(), categories=["prior malignancy"])
```

## Which entities are checked

`/annotations` matches an entity's **own** UUID: an annotation on a case is not
returned for that case's aliquots, nor the other way round. So
{meth}`~gdc2maf.cohort.Cohort.annotations` queries all three levels the cohort
depends on — the selected MAF files, their tumour and normal aliquots, and the
included cases — and saves the lot as `<name>_gdc_annotations.tsv`.

On the three documented cohorts, no selected file and no selected aliquot
carried an exclusion-grade annotation. The step is still worth running: it is
the record that the check was made, and a later GDC release can change the
answer.

## The PanCanAtlas supplements

{mod}`gdc2maf.pancan` downloads the two PanCanAtlas publication supplements and
verifies them:

| Name | File | Holds |
|---|---|---|
| `clinical` | `clinical_PANCAN_patient_with_followup.tsv` | 746 columns of curated clinical and follow-up data per patient |
| `quality_annotations` | `merged_sample_quality_annotations.tsv` | the `Do_not_use` flag and the AWG pathology exclusions |

```python
from gdc2maf import fetch_pancan_clinical, fetch_pancan_quality_annotations

record = fetch_pancan_clinical("data/TCGA")
cohort = Cohort(..., sex_fallback_path=record["path"])
```

These are publication supplements rather than indexed GDC data files — the
`/files` endpoint returns nothing for their UUIDs — so the API publishes no
checksum for them. {data}`gdc2maf.pancan.PANCAN_FILES` pins the md5 of the file
served, and it is checked on every download, so a file that changes upstream
fails loudly instead of quietly changing an analysis. A copy already on disk
with the right md5 is reused, so the call is cheap to leave in a pipeline.

For the clinical table, the GDC's own {func}`gdc2maf.clinical.fetch_clinical` is
the better default: it is pinned to the release the cohort was built from, and
its 33 columns already carry the usual variables. Reach for the PanCanAtlas one
when you need a column the GDC does not have.

## Reading the flag yourself

The pieces are usable without a `Cohort`:

```python
from gdc2maf import load_quality_annotations, do_not_use_patients

annotations = load_quality_annotations("data/TCGA/merged_sample_quality_annotations.tsv")
flagged = do_not_use_patients(annotations)                  # any platform
exome_only = do_not_use_patients(annotations, platform="WXS")
```

`do_not_use_patients` takes the conservative reading the TCGA analysis working
groups applied — a patient flagged on *any* platform is flagged — because the
flag is recorded per platform and a reason such as prior treatment applies to
the patient however the DNA was sequenced. Narrow it with `platform` to flag on
exome QC alone.
