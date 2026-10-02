# Tutorials

Worked examples, from the shortest useful run to a complete analysis pipeline.

:::{list-table}
:header-rows: 1
:widths: 35 65

* - Tutorial
  - What it covers
* - {doc}`first-cohort`
  - One cohort from scratch: find the projects, resolve the cases, download and
    verify the files, merge, and read the attrition report.
* - {doc}`signature-pipeline`
  - A mutational-signature pipeline built on the API: cohort to merged MAF with
    gdc2maf, then SigProfilerMatrixGenerator and SigProfilerExtractor, with
    per-cohort logging, resumable steps and the patients lost at each step
    accounted for.
:::

```{toctree}
:hidden:

first-cohort
signature-pipeline
```
