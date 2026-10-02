# Cohorts and file specs

A cohort is defined by GDC project(s) plus optional case-level criteria. This page
covers how to find the projects, what the criteria do, and which GDC files count as
the cohort's MAFs.

(finding-the-project-ids)=
## Finding the project IDs

`--project` takes GDC project IDs such as `TCGA-LUAD` — not a tissue or a disease name.
Two ways to look one up, both repeated in `gdc2maf -h`:

1. Browse <https://portal.gdc.cancer.gov/projects> and filter by primary site or program;
   the **Project** column holds the IDs.
2. Ask the API. Every project holding a lung case, as a table:

   ```bash
   curl -G https://api.gdc.cancer.gov/projects \
     --data-urlencode 'fields=project_id,name' \
     --data-urlencode 'size=100' \
     --data-urlencode 'format=TSV' \
     --data-urlencode 'filters={"op":"in","content":{"field":"primary_site","value":["Bronchus and lung"]}}'
   ```

   Drop the last line to list every project. Site names follow ICD-O-3, so lung is
   `Bronchus and lung`, not `Lung`.

**A primary site is not a disease.** That lung query returns 28 projects — `TCGA-MESO`
(mesothelioma) and `TCGA-SKCM` (melanoma) among them — because each holds at least one case
recorded at a lung site. This is why gdc2maf asks for projects rather than a site: pick the
projects for the disease you mean and pass them all to `--project`. TCGA lung cancer
is `TCGA-LUAD` (adenocarcinoma) plus `TCGA-LUSC` (squamous cell):

```bash
gdc2maf maf --project TCGA-LUAD TCGA-LUSC --name Lung
```

Every case of the named projects is fetched and kept in the cases table, so the attrition
report starts from the whole project and shows what each criterion removed.

(sex-criterion)=
## Choosing which cases: sex at birth

The one case-level criterion gdc2maf applies is `demographic.sex_at_birth`. It is
optional, takes one or several values, and is spelled the same way in both interfaces:

::::{tab-set}
:::{tab-item} Command line
```bash
gdc2maf maf --project TCGA-BRCA --sex-at-birth female --name Breast_FEMALE
gdc2maf maf --project TCGA-BRCA --sex-at-birth male   --name Breast_MALE
gdc2maf maf --project TCGA-BRCA --sex-at-birth male female --name Breast_BOTH
gdc2maf maf --project TCGA-BRCA --name Breast_ALL      # no criterion: every case
```
Several values after one flag, or the flag repeated; the two forms add up.
:::
:::{tab-item} Python
```python
Cohort(..., sex_at_birth=["female"])          # female only
Cohort(..., sex_at_birth=["male"])            # male only
Cohort(..., sex_at_birth=["male", "female"])  # both, explicitly
Cohort(..., sex_at_birth=None)                # no criterion: every case (the default)
```
:::
::::

### What each choice keeps

:::{list-table}
:header-rows: 1
:widths: 40 20 20 20

* - Criterion
  - male kept
  - female kept
  - sex not recorded
* - `["female"]` / `--sex-at-birth female`
  - no
  - yes
  - **no**
* - `["male"]` / `--sex-at-birth male`
  - yes
  - no
  - **no**
* - `["male", "female"]` / `--sex-at-birth male female`
  - yes
  - yes
  - **no**
* - `None` / the flag omitted
  - yes
  - yes
  - **yes**
:::

The last two rows are the ones to be deliberate about. Naming both sexes explicitly still
**excludes cases whose sex the GDC does not record** — they match neither value — whereas
leaving the criterion out applies no sex filter at all and keeps them. On a real cohort
that is not a rounding error: a male TCGA lung cohort lost 63 cases to a missing sex
against 411 to being female.

Nothing disappears quietly either way. Excluded cases stay in the cases table, marked with
the step and the reason — `sex_at_birth is female`, or
`sex_at_birth missing in GDC` — and are counted in the attrition report. See
{func}`gdc2maf.cases.fetch_cases`.

### Filling a missing sex from another table

Where the GDC records no sex for a case but another clinical table does, give that table
with `--sex-fallback` (Python: `sex_fallback_path=`) and it is used to fill the gap
*before* the criterion is applied, so those cases can still qualify. Only missing values
are filled; where both sources have a value the GDC's is kept and disagreements are
counted and logged. The defaults match the TCGA PanCan clinical table; see
{func}`gdc2maf.cases.fill_sex_from_table`.

## Which files are "the cohort's MAFs"

By default, the open-access GRCh38 `Masked Somatic Mutation` files from the GDC's Aliquot
Ensemble Somatic Variant Merging and Masking workflow, exome (WXS) — one MAF per
tumour-normal aliquot pair. That is a value, not hard-coded logic:

```python
from gdc2maf import Cohort, FileSpec

wgs_ensemble = FileSpec(
    name="wgs_ensemble_maf",
    title="open WGS ensemble MAF",
    data_category=("Simple Nucleotide Variation",),
    data_type=("Masked Somatic Mutation",),
    experimental_strategy=("WGS",),
    access=("open",),
)
cohort = Cohort(..., spec=wgs_ensemble)
```

Only `WXS_ENSEMBLE_MAF` is exercised by this package's own runs; other specs are passed
straight through to the GDC `/files` endpoint.

