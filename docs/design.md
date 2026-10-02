# Design notes


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
