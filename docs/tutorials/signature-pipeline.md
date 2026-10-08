# A mutational-signature pipeline

gdc2maf stops at the merged MAF. This adds the two
[SigProfiler](https://github.com/AlexandrovLab) steps on top, so one script takes
you from a cohort definition to extracted mutational signatures.

The cohort here is male lung-cancer patients in TCGA, which is two projects —
`TCGA-LUAD` and `TCGA-LUSC`. {ref}`finding-the-project-ids` shows how to look up
the IDs for a different disease.

| Stage | Done by |
|---|---|
| Cohort → verified, merged MAF | gdc2maf, in one call |
| Merged MAF → SigProfiler input | a few lines of reshaping |
| Input → mutation matrices | SigProfilerMatrixGenerator |
| Matrices → signatures | SigProfilerExtractor |

## Before you start

```bash
pip install gdc2maf SigProfilerMatrixGenerator SigProfilerExtractor
```

SigProfiler needs its reference genome installed once, which is a few GB:

```python
from SigProfilerMatrixGenerator import install as genInstall

genInstall.install("GRCh38")
```

## The pipeline

```python
"""Mutational signatures for male lung-cancer patients in TCGA."""

import os

from SigProfilerExtractor import sigpro
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

from gdc2maf import Cohort, configure_logging

# SigProfilerMatrixGenerator reads MAF columns by position, so the input must
# keep the first 16 GDC columns in their original order.
MAF_COLUMNS = [
    "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome",
    "Start_Position", "End_Position", "Strand", "Variant_Classification",
    "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
    "dbSNP_RS", "dbSNP_Val_Status", "Tumor_Sample_Barcode",
]
# It treats anything longer than two bases as an indel, so multi-base
# substitutions (TNP, ONP) would land in the indel matrix. Keep only these.
VARIANT_TYPES = ["SNP", "DNP", "INS", "DEL"]


def write_sigprofiler_input(cohort, sigprofiler_dir):
    """Write the MatrixGenerator input MAF from the cohort's merged MAF."""
    maf = cohort.read_merged_maf(usecols=MAF_COLUMNS)
    maf = maf[maf["Variant_Type"].isin(VARIANT_TYPES)]

    input_dir = f"{sigprofiler_dir}/input"
    os.makedirs(input_dir, exist_ok=True)
    maf.to_csv(f"{input_dir}/{cohort.name}.maf", sep="\t", index=False)
    print(f"SigProfiler input: {len(maf)} variants, "
          f"{maf['Tumor_Sample_Barcode'].nunique()} samples")


def main():
    """Run the whole pipeline: cohort, MAF, matrices, signatures."""
    configure_logging()

    # Cohort -> one verified, merged MAF, with QC and an attrition report
    cohort = Cohort(
        name="Lung_MALE",
        projects=["TCGA-LUAD", "TCGA-LUSC"],
        sex_at_birth=["male"],
        out_dir="out/Lung_MALE",
        download_dir="downloads",   # shared by every cohort; see ../outputs
    )
    cohort.run()

    sigprofiler_dir = f"{cohort.out_dir}/sigprofiler"
    write_sigprofiler_input(cohort, sigprofiler_dir)

    # Mutation matrices (minutes), written to {sigprofiler_dir}/output/
    matGen.SigProfilerMatrixGeneratorFunc(
        cohort.name, "GRCh38", sigprofiler_dir, exome=True
    )

    # Signature extraction (hours). For matrix input the Extractor uses
    # opportunity_genome; reference_genome only applies to vcf input.
    sigpro.sigProfilerExtractor(
        "matrix",
        f"{sigprofiler_dir}/extractor",
        f"{sigprofiler_dir}/output/SBS/{cohort.name}.SBS96.exome",
        opportunity_genome="GRCh38",
        exome=True,
        cosmic_version=3.6,
        minimum_signatures=1,
        maximum_signatures=20,
        nmf_replicates=100,
    )


if __name__ == "__main__":
    main()
```

## What to watch out for

**Column order is not cosmetic.** MatrixGenerator's MAF reader takes columns by
*position*, not by name: Chromosome 5, Start_Position 6, End_Position 7,
Reference_Allele 11, Tumor_Seq_Allele2 13, Tumor_Sample_Barcode 16. `MAF_COLUMNS`
above is the GDC MAF's own first 16 columns, which line up. Reorder them and the
matrices come out silently wrong rather than empty.

**`input/` should hold nothing else.** MatrixGenerator reads every file in that
folder, so if you rerun with a different cohort name, clear it first.

**The merged MAF has already had the flagged patients removed.** `cohort.run()`
leaves out the patients the GDC's curation or the PanCanAtlas `Do_not_use` flag
names — 56 of 596 for this cohort — so the matrices and the extracted signatures
never see them, and `Lung_MALE_case_attrition.tsv` reports them at step 6.
`exclude_flagged=False` keeps them. Either way the choice belongs in the
write-up of the analysis, so read {doc}`../quality` once before running the
extraction; changing it afterwards means re-merging (`refresh=["merge"]`) and
rerunning SigProfiler on new matrices.

**One `download_dir` for all your cohorts.** A study is rarely one cohort, and
`downloads/` is keyed by GDC file UUID, not by cohort, so the same directory
serves every one of them and a file is fetched once. Running this script for a
second cohort with its own download directory would fetch that cohort's files
again from scratch whenever the two overlap — a sensitivity subset is the usual
case. {doc}`../outputs` has the detail.

**`read_merged_maf(usecols=...)` reads only what you need.** A cohort's merged MAF
runs to hundreds of MB; this pulls back the 16 columns instead of all of them.

## Where the outputs land

```text
out/Lung_MALE/
├── Lung_MALE_cases.tsv                     gdc2maf
├── Lung_MALE_files_selected.tsv
├── Lung_MALE_download_check.tsv
├── Lung_MALE_gdc_annotations.tsv
├── Lung_MALE_quality_flags.tsv
├── Lung_MALE_wxs_ensemble_maf.maf
├── Lung_MALE_sample_qc.tsv
├── Lung_MALE_case_attrition.tsv
└── sigprofiler/                            SigProfiler
    ├── input/Lung_MALE.maf
    ├── output/SBS/Lung_MALE.SBS96.exome
    └── extractor/SBS96/Suggested_Solution/
```

{doc}`../outputs` describes the gdc2maf files; the SigProfiler layout is
documented upstream.

## Running it again

Rerunning the script does not redo the work: every gdc2maf step reuses what it
left in `out_dir`, and files already downloaded are verified rather than fetched
again. SigProfiler is the slow half, and it does not cache, so once the matrices
exist drop the `SigProfilerMatrixGeneratorFunc(...)` call from `main()`, and once
a solution exists drop the `sigProfilerExtractor(...)` call as well.

For a long extraction a separate job is usually better than one long script. The
matrices are already on disk by then, so that job is just the last call:

```python
from SigProfilerExtractor import sigpro

sigprofiler_dir = "out/Lung_MALE/sigprofiler"
sigpro.sigProfilerExtractor(
    "matrix",
    f"{sigprofiler_dir}/extractor",
    f"{sigprofiler_dir}/output/SBS/Lung_MALE.SBS96.exome",
    opportunity_genome="GRCh38",
    exome=True,
    cosmic_version=3.6,
    minimum_signatures=1,
    maximum_signatures=20,
    nmf_replicates=100,
)
```

{doc}`../outputs` explains which gdc2maf steps are cached and how `refresh` asks
for one to be redone; {doc}`../cli` has the per-step subcommands for splitting the
gdc2maf half the same way.

## Reporting the patients SigProfiler drops

Step 2 keeps only four variant types, so a patient whose variants are all TNP or
ONP has variants in the merged MAF but none in the SigProfiler input. To account
for those in the same attrition table as the GDC steps, pass them to
{func}`~gdc2maf.reports.summarize_attrition` as an extra step — it takes
`(label, case_ids_lost, reason)` triples, so a step gdc2maf knows nothing about
still appears in the report.

gdc2maf's own steps are numbered 1 to 7, so a downstream step starts at 8:

```python
STEP_SIGPROFILER_INPUT = "8. SigProfiler input"

with_variants = sample_qc[sample_qc["has_variants"].astype(bool)]
not_in_input = with_variants.loc[
    ~with_variants["Tumor_Sample_Barcode"].isin(
        sigprofiler_input["Tumor_Sample_Barcode"]
    ),
    "case_id",
]
attrition = cohort.attrition(
    extra_steps=[
        (STEP_SIGPROFILER_INPUT, not_in_input, "no variants of the kept types")
    ]
)
```
