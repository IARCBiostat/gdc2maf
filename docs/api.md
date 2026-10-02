# Python API

[The quick start](index.md#quick-start) shows the whole run in one call,
{meth}`cohort.run() <gdc2maf.cohort.Cohort.run>`. This page does the same work step by
step. Each method returns DataFrames and caches to `out_dir`; calling a late step runs the
earlier ones it needs, so you can start anywhere.

```{important}
`projects` takes GDC **project IDs** such as `TCGA-LUAD`, not a tissue or disease name,
and a disease usually spans more than one project. Start at
{ref}`finding-the-project-ids` to pick the right ones from the GDC portal or its API —
filtering by tissue alone sweeps in projects that merely contain a case biopsied there.
```

```python
from gdc2maf import Cohort, configure_logging

configure_logging()  # the package logs nothing until you ask

cohort = Cohort(
    name="Lung_MALE",
    projects=["TCGA-LUAD", "TCGA-LUSC"],
    sex_at_birth=["male"],
    out_dir="out/Lung_MALE",
    download_dir="downloads",
)

cases = cohort.cases()  # GDC query, or the cached TSV
maf_files = cohort.maf_files()  # one row per candidate file, nothing downloaded
selected, dropped = cohort.selection()  # one file per patient + why the others went
check = cohort.download(n_clients=8)  # gdc-client, then size + md5 on every file
maf, sample_qc, qc = cohort.merged_maf()

print(cohort.attrition().to_string(index=False))
```

{meth}`cohort.run() <gdc2maf.cohort.Cohort.run>` does all of it and returns every table in
a dict — that is the form [the quick start](index.md#quick-start) uses. To work from a
merged MAF that is already on disk without loading all of it, use
{meth}`cohort.read_merged_maf(usecols=[...]) <gdc2maf.cohort.Cohort.read_merged_maf>`.

## Choosing the cases

`sex_at_birth` takes a list of values, or `None` to apply no criterion:

```python
Cohort(..., sex_at_birth=["female"])          # female only
Cohort(..., sex_at_birth=["male"])            # male only
Cohort(..., sex_at_birth=["male", "female"])  # both, explicitly
Cohort(..., sex_at_birth=None)                # no criterion: every case (the default)
```

The last two are not the same — naming both sexes still excludes cases whose sex the GDC
does not record, while `None` keeps them. What each choice keeps, and how
`sex_fallback_path` fills a missing sex, is set out once in {ref}`sex-criterion`.

The step functions are importable on their own (`fetch_cases`, `list_maf_files`,
`select_one_file_per_case`, `download_files`, `merge_mafs`, `ensure_gdc_client`, …) if you
would rather drive them yourself than hold a `Cohort`.


See {doc}`reference/index` for every class and function, and
{doc}`tutorials/index` for worked examples.
