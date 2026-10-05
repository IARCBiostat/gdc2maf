# gdc2maf

```{include} ../README.md
:parser: myst_parser.sphinx_
:start-after: "<!-- docs-include: intro-start -->"
:end-before: "<!-- docs-include: intro-end -->"
```

```{include} ../README.md
:parser: myst_parser.sphinx_
:start-after: "<!-- docs-include: quickstart-start -->"
:end-before: "<!-- docs-include: quickstart-end -->"
```

## Where to go next

:::{list-table}
:header-rows: 1
:widths: 30 70

* - Page
  - What is on it
* - {doc}`install`
  - Installing with pip or conda, and which platforms `gdc-client` supports.
* - {doc}`cli`
  - Every command and option, how to split a run across cluster jobs, and how to
    find the GDC project IDs for the disease you mean.
* - {doc}`api`
  - The same steps driven from Python, one at a time, with what each returns.
* - {doc}`cohorts`
  - What defines a cohort, and which GDC files count as its MAFs.
* - {doc}`outputs`
  - Every file a run writes, and exactly what a rerun does and does not redo.
* - {doc}`quality`
  - The GDC's curation annotations and the PanCanAtlas `Do_not_use` flag: what
    is checked, what is excluded by default, and when to keep it instead.
* - {doc}`tutorials/index`
  - Worked examples: a first cohort end to end, and a full mutational-signature
    pipeline built on the API.
* - {doc}`reference/index`
  - Every class and function, generated from the docstrings.
* - {doc}`design`
  - Why the package behaves as it does: what it refuses, what it only flags.
:::

## The API at a glance

Everything marked below is importable straight from `gdc2maf`
(`from gdc2maf import Cohort`). Follow a link for the full signature and
parameters, or see {doc}`reference/index` for every object grouped by module.

### The usual way in

```{eval-rst}
.. autosummary::
   :nosignatures:

   gdc2maf.cohort.Cohort
   gdc2maf.provenance.configure_logging
```

### Which files are the cohort's MAFs

```{eval-rst}
.. autosummary::
   :nosignatures:

   gdc2maf.spec.FileSpec
   gdc2maf.spec.WXS_ENSEMBLE_MAF
```

### The steps, one at a time

```{eval-rst}
.. autosummary::
   :nosignatures:

   gdc2maf.cases.fetch_cases
   gdc2maf.clinical.fetch_clinical
   gdc2maf.files.list_maf_files
   gdc2maf.files.inspect_duplicates
   gdc2maf.selection.select_one_file_per_case
   gdc2maf.client.ensure_gdc_client
   gdc2maf.download.download_files
   gdc2maf.maf.merge_mafs
```

### Reports, records and logging

```{eval-rst}
.. autosummary::
   :nosignatures:

   gdc2maf.reports.summarize_attrition
   gdc2maf.reports.summary_table
   gdc2maf.record.download_record_text
   gdc2maf.provenance.log_environment
   gdc2maf.provenance.tee_stdout
```

```{toctree}
:hidden:
:caption: Using gdc2maf

install
cli
api
cohorts
outputs
quality
```

```{toctree}
:hidden:
:caption: Examples

tutorials/index
```

```{toctree}
:hidden:
:caption: Reference

reference/index
design
```
