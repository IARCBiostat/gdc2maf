# Command line

## One command does the whole job

`gdc2maf maf` runs every step in order — resolve the cohort, list its MAF files, pick one
per patient, install `gdc-client`, download and verify, merge, and write the reports. For
most work it is the only command you need:

```bash
gdc2maf maf --project TCGA-LUAD TCGA-LUSC \
            --sex-at-birth male --name Lung_MALE \
            --out out --download-dir downloads --jobs 8
```

That writes everything to `out/Lung_MALE/`, downloads into `downloads/<file_id>/`, and
prints the attrition table and a summary of each step.

```{important}
`--project` takes GDC **project IDs** such as `TCGA-LUAD`, not a tissue or disease name,
and a disease usually spans more than one project. Start at
{ref}`finding-the-project-ids` to pick the right ones from the GDC portal or its API —
filtering by tissue alone sweeps in projects that merely contain a case biopsied there.
```

Each step is *also* available as its own subcommand, for splitting a long run across
cluster jobs or redoing one stage on its own — see
[the step subcommands](#the-step-subcommands) at the end of this page.

## Choosing the cases

`--sex-at-birth` restricts the cohort by `demographic.sex_at_birth`. Several values after
one flag, or the flag repeated:

```bash
gdc2maf maf --project TCGA-BRCA --sex-at-birth female --name Breast_FEMALE
gdc2maf maf --project TCGA-BRCA --sex-at-birth male   --name Breast_MALE
gdc2maf maf --project TCGA-BRCA --sex-at-birth male female --name Breast_BOTH
gdc2maf maf --project TCGA-BRCA --name Breast_ALL      # no criterion: every case
```

The last two are not the same — naming both sexes still excludes cases whose sex the GDC
does not record, while omitting the flag keeps them. What each choice keeps, and how
`--sex-fallback` fills a missing sex, is set out once in
{ref}`sex-criterion`.

## Naming projects

`--project` takes as many IDs as you like after one flag, and the flag can also be
repeated; the two forms add up, so these are all the same cohort:

```bash
--project TCGA-LUAD TCGA-LUSC
--project TCGA-LUAD --project TCGA-LUSC
--project TCGA-LUAD TCGA-LUSC --project TCGA-KIRC   # three projects
```

The same goes for the other repeatable options — `--sex-at-birth`, `--refresh` and
`--remove-gdc-filter-flag`. For how to find the IDs for the disease you mean, see
{doc}`cohorts`.

## The gdc-client is handled for you

Downloading GDC files needs the GDC Data Transfer Tool, so before the first download gdc2maf looks for `gdc-client` in
`tools/gdc-client/<version>/` and, if it is not there, downloads that release from the GDC
and checks its md5 against a pinned table before running it. Which binary ran — path,
version, and its own md5 — is recorded in `<name>_gdc_client.json` beside the cohort's
other outputs. Version `2.3` is pinned by default; `--gdc-client-version latest` takes the
newest release, `--gdc-client PATH` uses a client you already have, and `--install-dir`
moves where versions are kept. `gdc2maf client` does only this step, which is useful for
warming a shared installation before a batch of cluster jobs.

## Rerunning

A rerun reuses what is already on disk, so it is cheap to run again after an interruption:

```bash
gdc2maf maf ... --refresh files      # re-query the file listing, keep the rest
gdc2maf maf ... --no-download        # verify the files on disk, download nothing
```

{doc}`outputs` has the full account of what a rerun does and does not redo.

## The step subcommands

Everything above runs inside `gdc2maf maf`. When you want the stages separately — a
cluster job per step, or to redo one without the others — each is its own subcommand,
taking the same cohort options:

| Command | Does |
|---|---|
| `gdc2maf cases` | Resolve the cohort into a table of GDC cases |
| `gdc2maf clinical` | Fetch clinical data for those cases (auxiliary) |
| `gdc2maf files` | List the MAF files; describe patients with more than one |
| `gdc2maf select` | Choose one file per patient, with the reason each other was dropped |
| `gdc2maf download` | Download the selected files; verify size and md5 |
| `gdc2maf merge` | Merge the verified MAFs into one MAF, with QC |
| `gdc2maf reports` | Rebuild the attrition and summary tables from cached outputs |
| `gdc2maf client` | Locate or install `gdc-client` and report which binary is used |
| `gdc2maf maf` | All of the above, in order |

A step reuses whatever the previous run left in the output folder, so running them one
after another is equivalent to the single `gdc2maf maf` call. See
{doc}`outputs` for what each step writes and what a rerun redoes.

## Help and common options

Every command has built-in help: `gdc2maf -h` (or `--help`) lists the commands, and
`gdc2maf <command> -h` lists that command's options with their defaults — reproduced in
full under [Option reference](#option-reference) below. Useful ones: `--sex-fallback`
(fill `sex_at_birth` from an external clinical table where the GDC has none),
`--remove-gdc-filter-flag` (drop `GDC_FILTER`-flagged variants; the default keeps them and
reports the counts), `--gdc-client-version` (pin a version, or `latest`), `--no-log`,
`--quiet`.


For how to find the project IDs to pass, and what the cohort criteria do, see
{doc}`cohorts`.

## Option reference

Generated from the parser itself, so it cannot drift from what `gdc2maf -h` prints.

```{eval-rst}
.. argparse::
   :module: gdc2maf.cli
   :func: build_parser
   :prog: gdc2maf
```
