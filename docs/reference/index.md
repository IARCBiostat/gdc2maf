# API reference

Generated from the docstrings, so nothing here is written twice: the source is
the only place a signature or a parameter is described.

Every object below is documented on its module's page. The names marked
*(top level)* are also importable straight from `gdc2maf`, which is the usual
way to reach them — `from gdc2maf import Cohort`.

## Modules

```{eval-rst}
.. autosummary::
   :toctree: generated
   :recursive:

   gdc2maf.cohort
   gdc2maf.spec
   gdc2maf.cases
   gdc2maf.clinical
   gdc2maf.files
   gdc2maf.selection
   gdc2maf.client
   gdc2maf.download
   gdc2maf.maf
   gdc2maf.attrition
   gdc2maf.annotations
   gdc2maf.pancan
   gdc2maf.reports
   gdc2maf.record
   gdc2maf.provenance
   gdc2maf.tcga
   gdc2maf.api
   gdc2maf.cli
```

## By topic

### The cohort object *(top level)*

```{eval-rst}
.. autosummary::
   :nosignatures:

   gdc2maf.cohort.Cohort
   gdc2maf.cohort.CACHED_STEPS
```

### Which files are the cohort's MAFs *(top level)*

```{eval-rst}
.. autosummary::
   :nosignatures:

   gdc2maf.spec.FileSpec
   gdc2maf.spec.WXS_ENSEMBLE_MAF
```

### Resolving a cohort *(top level)*

```{eval-rst}
.. autosummary::
   :nosignatures:

   gdc2maf.cases.fetch_cases
   gdc2maf.cases.fill_sex_from_table
   gdc2maf.clinical.fetch_clinical
   gdc2maf.files.list_maf_files
   gdc2maf.files.inspect_duplicates
   gdc2maf.selection.select_one_file_per_case
   gdc2maf.tcga.parse_tcga_aliquot_barcode
```

### Downloading *(top level)*

```{eval-rst}
.. autosummary::
   :nosignatures:

   gdc2maf.client.ensure_gdc_client
   gdc2maf.client.install_gdc_client
   gdc2maf.download.write_manifest
   gdc2maf.download.check_downloaded_files
   gdc2maf.download.download_files
   gdc2maf.download.file_md5
```

### Merging and QC

```{eval-rst}
.. autosummary::
   :nosignatures:

   gdc2maf.maf.merge_mafs
   gdc2maf.maf.read_gdc_maf
   gdc2maf.maf.sample_qc_metrics
   gdc2maf.maf.flag_samples
```

### Attrition and reports

```{eval-rst}
.. autosummary::
   :nosignatures:

   gdc2maf.attrition.exclude_cases
   gdc2maf.attrition.maf_step
   gdc2maf.reports.summarize_attrition
   gdc2maf.reports.summarize_maf_files
   gdc2maf.reports.summarize_duplicates
   gdc2maf.reports.summarize_selection
   gdc2maf.reports.summarize_download
   gdc2maf.reports.summary_table
```

### Quality flags: GDC curation and the PanCanAtlas

```{eval-rst}
.. autosummary::
   :nosignatures:

   gdc2maf.annotations.fetch_annotations
   gdc2maf.annotations.exclusion_flags
   gdc2maf.annotations.summarize_annotations
   gdc2maf.annotations.EXCLUSION_CLASSIFICATIONS
   gdc2maf.annotations.EXCLUSION_CATEGORIES
   gdc2maf.pancan.fetch_pancan_file
   gdc2maf.pancan.fetch_pancan_clinical
   gdc2maf.pancan.fetch_pancan_quality_annotations
   gdc2maf.pancan.load_quality_annotations
   gdc2maf.pancan.do_not_use_patients
   gdc2maf.pancan.write_pancan_record
   gdc2maf.pancan.PANCAN_FILES
```

### The download record

```{eval-rst}
.. autosummary::
   :nosignatures:

   gdc2maf.record.download_record_text
   gdc2maf.record.write_if_changed
```

### Logging and run records *(top level)*

```{eval-rst}
.. autosummary::
   :nosignatures:

   gdc2maf.provenance.configure_logging
   gdc2maf.provenance.log_environment
   gdc2maf.provenance.environment_info
   gdc2maf.provenance.tee_stdout
```

### The GDC API layer

```{eval-rst}
.. autosummary::
   :nosignatures:

   gdc2maf.api.gdc_query
   gdc2maf.api.gdc_data_release
   gdc2maf.api.gdc_in
   gdc2maf.api.gdc_and
   gdc2maf.api.gdc_cohort_cases
   gdc2maf.api.cases_with_files
   gdc2maf.api.maf_file_metadata
   gdc2maf.clinical.gdc_clinical
   gdc2maf.clinical.gdc_case_country
```

### Command line internals

```{eval-rst}
.. autosummary::
   :nosignatures:

   gdc2maf.cli.build_parser
   gdc2maf.cli.main
   gdc2maf.cli.cohort_from_args
   gdc2maf.cli.default_cohort_name
```
