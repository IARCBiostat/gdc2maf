"""Command line interface: ``gdc2maf``.

``gdc2maf maf`` runs the whole path from a cohort specification to a verified,
merged MAF. The individual steps are also subcommands, so a long run can be
split across jobs (a cluster job per cohort, or per step) and reruns reuse what
is already on disk.
"""

import argparse
import logging
import sys

from . import __version__
from .client import ensure_gdc_client
from .cohort import CACHED_STEPS, Cohort
from .provenance import configure_logging, log_environment
from .reports import summary_table

logger = logging.getLogger(__name__)

GDC_PROJECTS_PORTAL = "https://portal.gdc.cancer.gov/projects"
GDC_PROJECTS_API = "https://api.gdc.cancer.gov/projects"

# One long line on purpose: it has to survive copy-and-paste out of the help.
_SITE_FILTER = (
    "--data-urlencode 'filters={\"op\":\"in\",\"content\":"
    "{\"field\":\"primary_site\",\"value\":[\"Bronchus and lung\"]}}'"
)

#: Printed after the options of every command that takes ``--project``.
FINDING_PROJECTS = (
    "finding project IDs:\n"
    "  --project takes GDC project IDs such as TCGA-LUAD, not a tissue or a\n"
    "  disease name. Two ways to look one up:\n"
    "\n"
    f"  1. Browse {GDC_PROJECTS_PORTAL} and filter by primary\n"
    '     site or program; the "Project" column holds the IDs.\n'
    "\n"
    "  2. Ask the API. Every project holding a lung case, as a table:\n"
    "\n"
    f"       curl -G {GDC_PROJECTS_API} \\\n"
    "         --data-urlencode 'fields=project_id,name' \\\n"
    "         --data-urlencode 'size=100' \\\n"
    "         --data-urlencode 'format=TSV' \\\n"
    f"         {_SITE_FILTER}\n"
    "\n"
    "     Drop the last line to list every project. Site names follow\n"
    '     ICD-O-3, so lung is "Bronchus and lung", not "Lung".\n'
    "\n"
    "  Mind that a primary site is not a disease: that lung query returns 28\n"
    "  projects, TCGA-MESO (mesothelioma) and TCGA-SKCM (melanoma) among them,\n"
    "  because each holds some case recorded at a lung site. Pick the projects\n"
    "  for the disease you mean and pass them all to --project. TCGA lung\n"
    "  cancer, for instance, is TCGA-LUAD (adenocarcinoma) and TCGA-LUSC\n"
    "  (squamous cell):\n"
    "\n"
    "       gdc2maf maf --project TCGA-LUAD TCGA-LUSC --name Lung\n"
)


def default_cohort_name(projects, sex_at_birth):
    """Build a cohort name from the projects and sex, for use without ``--name``."""
    name = "_".join(projects)
    if sex_at_birth:
        name += "_" + "_".join(s.upper() for s in sex_at_birth)
    return name


def add_cohort_arguments(parser):
    """Add the arguments that define a cohort and where its outputs go."""
    group = parser.add_argument_group("cohort")
    group.add_argument(
        "--project",
        action="extend",
        nargs="+",
        required=True,
        metavar="ID",
        help="GDC project ID, e.g. TCGA-LUAD. Several at once (--project "
        "TCGA-LUAD TCGA-LUSC) or the flag repeated (--project TCGA-LUAD "
        "--project TCGA-LUSC); both add up. Every case of these projects is "
        "fetched and reported, so the cohort's attrition starts from the whole "
        "project. See 'finding project IDs' at the end of this help for how to "
        f"look an ID up ({GDC_PROJECTS_PORTAL})",
    )
    group.add_argument(
        "--sex-at-birth",
        action="extend",
        nargs="+",
        choices=["male", "female"],
        metavar="SEX",
        help="keep only cases with this demographic.sex_at_birth; several at "
        "once (--sex-at-birth male female) or the flag repeated. Choices: male, "
        "female. Default: every sex is kept, including cases whose sex the GDC "
        "does not record",
    )
    group.add_argument(
        "--name",
        help="cohort label used in every output file name "
        "(default: the projects and sex joined)",
    )
    group.add_argument(
        "--sex-fallback",
        metavar="TSV",
        help="clinical table used to fill sex_at_birth where the GDC has none "
        "(e.g. the TCGA PanCan clinical table)",
    )
    group.add_argument(
        "--sex-fallback-source",
        default="PanCan",
        metavar="LABEL",
        help="label recorded for cases filled from --sex-fallback "
        "(default: %(default)s)",
    )

    group = parser.add_argument_group("quality flags")
    group.add_argument(
        "--quality-annotations",
        metavar="TSV",
        help="an existing PanCanAtlas merged_sample_quality_annotations.tsv to "
        "read the Do_not_use flag from. Omitted, the table is downloaded "
        "(md5-checked) into --pancan-dir the first time it is needed",
    )
    group.add_argument(
        "--pancan-dir",
        default="data/pancan",
        metavar="DIR",
        help="where the PanCanAtlas table is downloaded and reused from; safe "
        "to share between cohorts (default: %(default)s)",
    )
    group.add_argument(
        "--no-quality-annotations",
        dest="fetch_quality_annotations",
        action="store_false",
        help="do not download or read the PanCanAtlas table; only the GDC's "
        "own curation annotations can then flag a patient",
    )
    group.add_argument(
        "--keep-flagged",
        dest="exclude_flagged",
        action="store_false",
        help="keep the patients the quality flags list instead of excluding "
        "them; <name>_quality_flags.tsv is written either way. Default: drop "
        "them and record the loss in the attrition table",
    )

    group = parser.add_argument_group("output")
    group.add_argument(
        "--out",
        default="out",
        metavar="DIR",
        help="parent output directory; this cohort writes to DIR/<name> "
        "(default: %(default)s)",
    )
    group.add_argument(
        "--download-dir",
        default="downloads",
        metavar="DIR",
        help="where gdc-client stores files, as DIR/<file_id>/<file_name>. Use "
        "ONE directory for every cohort: the folders are named after the GDC "
        "file UUID and the merge looks its files up by UUID, so a file is "
        "downloaded once however many cohorts select it. A directory per "
        "cohort downloads it again for each (default: %(default)s)",
    )
    group.add_argument(
        "--refresh",
        action="extend",
        nargs="+",
        choices=[*CACHED_STEPS, "all"],
        default=[],
        metavar="STEP",
        help="recompute this step instead of reusing its cached output; several "
        "at once (--refresh files download) or the flag repeated. Choices: "
        + ", ".join([*CACHED_STEPS, "all"]),
    )
    group.add_argument(
        "--no-log",
        action="store_true",
        help="do not write <name>_run.log (log to the terminal only)",
    )
    group.add_argument(
        "--quiet", action="store_true", help="log warnings and errors only"
    )


def add_download_arguments(parser):
    """Add the gdc-client arguments."""
    group = parser.add_argument_group("download")
    group.add_argument(
        "--jobs",
        type=int,
        default=8,
        metavar="N",
        help="parallel gdc-client processes (default: %(default)s)",
    )
    group.add_argument(
        "--no-download",
        action="store_true",
        help="verify the files already on disk and download nothing",
    )
    group.add_argument(
        "--gdc-client",
        metavar="PATH",
        help="use this gdc-client executable instead of installing one",
    )
    group.add_argument(
        "--gdc-client-version",
        default="2.3",
        metavar="V",
        help="gdc-client version to install and pin, or 'latest' "
        "(default: %(default)s)",
    )
    group.add_argument(
        "--record-dir",
        metavar="DIR",
        help="where to put <name>_download_record.txt, the plain-text record of "
        "how many cases were excluded and how many files were downloaded "
        "(default: next to the downloaded MAFs, in --download-dir)",
    )
    group.add_argument(
        "--install-dir",
        default="tools/gdc-client",
        metavar="DIR",
        help="where gdc-client versions are installed (default: %(default)s)",
    )


def add_merge_arguments(parser):
    """Add the merge arguments."""
    group = parser.add_argument_group("merge")
    group.add_argument(
        "--remove-gdc-filter-flag",
        action="extend",
        nargs="+",
        default=[],
        metavar="FLAG",
        help="drop variants carrying this GDC_FILTER flag; several at once "
        "(--remove-gdc-filter-flag NonExonic gdc_pon) or the flag repeated. "
        "Known flags: NonExonic, gdc_pon, common_in_gnomAD. Default: keep "
        "every flagged variant and report the counts",
    )


def build_parser():
    """Build the ``gdc2maf`` argument parser."""
    parser = argparse.ArgumentParser(
        prog="gdc2maf",
        description="Turn a GDC cohort specification into one verified, merged MAF.",
        epilog=FINDING_PROJECTS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"gdc2maf {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True, metavar="COMMAND")

    def sub(name, help_text, projects=True):
        """Add one subcommand; ``projects`` appends the project-ID guidance."""
        return subparsers.add_parser(
            name,
            help=help_text,
            description=help_text,
            epilog=FINDING_PROJECTS if projects else None,
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )

    p = sub(
        "maf", "Run the whole path: cases, files, selection, download, merge, reports."
    )
    add_cohort_arguments(p)
    add_download_arguments(p)
    add_merge_arguments(p)
    p.add_argument(
        "--no-clinical",
        action="store_true",
        help="skip the clinical table (it is not needed for the MAF)",
    )
    p.add_argument(
        "--no-merge",
        action="store_true",
        help="stop after the download check, without merging",
    )

    p = sub("cases", "Resolve the cohort into a table of GDC cases.")
    add_cohort_arguments(p)

    p = sub("clinical", "Fetch clinical data for the cohort's cases.")
    add_cohort_arguments(p)

    p = sub("files", "List the cohort's MAF files and inspect patients with several.")
    add_cohort_arguments(p)

    p = sub(
        "select",
        "Choose one MAF file per patient, with the reason each other was dropped.",
    )
    add_cohort_arguments(p)

    p = sub("download", "Download the selected files and verify size and md5.")
    add_cohort_arguments(p)
    add_download_arguments(p)

    p = sub("merge", "Merge the verified MAFs into one MAF, with QC.")
    add_cohort_arguments(p)
    add_merge_arguments(p)
    p.add_argument(
        "--record-dir",
        metavar="DIR",
        help="where to put <name>_download_record.txt "
        "(default: next to the downloaded MAFs, in --download-dir)",
    )

    p = sub("reports", "Rebuild the attrition and summary tables from cached outputs.")
    add_cohort_arguments(p)

    p = sub(
        "client",
        "Locate or install gdc-client and report which binary would be used.",
        projects=False,
    )
    add_download_arguments(p)

    return parser


def cohort_from_args(args):
    """Build a :class:`gdc2maf.cohort.Cohort` from parsed arguments."""
    name = args.name or default_cohort_name(args.project, args.sex_at_birth)
    return Cohort(
        name=name,
        projects=args.project,
        out_dir=f"{args.out}/{name}",
        download_dir=args.download_dir,
        sex_at_birth=args.sex_at_birth,
        sex_fallback_path=args.sex_fallback,
        sex_fallback_source=args.sex_fallback_source,
        remove_gdc_filter_flags=getattr(args, "remove_gdc_filter_flag", []),
        quality_annotations_path=args.quality_annotations,
        pancan_dir=args.pancan_dir,
        fetch_quality_annotations=args.fetch_quality_annotations,
        exclude_flagged=args.exclude_flagged,
        refresh=args.refresh,
    )


def log_table(title, table, index=True):
    """Log a DataFrame under a title, the way the reports are printed."""
    logger.info(f"\n{title}")
    logger.info(table.to_string(index=index))


def print_reports(cohort, extra_steps=()):
    """Log the cohort's attrition table and every available summary."""
    log_table(
        "Case attrition",
        cohort.attrition(extra_steps).drop(columns="cohort"),
        index=False,
    )
    titles = {
        "maf_files_summary": "MAF files",
        "duplicates_summary": (
            "Patients with >1 MAF file: how the files differ (flags overlap)"
        ),
        "selection_summary": "One file per patient",
        "download_summary": "Download check (size + md5)",
        "merge_qc_summary": "Merged MAF QC",
    }
    for report_name, summary in cohort.summaries().items():
        log_table(titles.get(report_name, report_name), summary_table([summary]))


def main(argv=None):
    """Entry point for the ``gdc2maf`` console script."""
    args = build_parser().parse_args(argv)

    if args.command == "client":
        configure_logging()
        info = ensure_gdc_client(
            version=args.gdc_client_version,
            install_dir=args.install_dir,
            gdc_client_path=args.gdc_client,
        )
        for key, value in info.items():
            logger.info(f"  {key}: {value}")
        return 0

    cohort = cohort_from_args(args)
    configure_logging(
        log_path=None if args.no_log else cohort.log_path,
        level=logging.WARNING if args.quiet else logging.INFO,
    )
    log_environment()

    if args.command == "maf":
        cohort.run(
            record_dir=args.record_dir,
            clinical=not args.no_clinical,
            download=not args.no_download,
            merge=not args.no_merge,
            gdc_client=args.gdc_client,
            n_clients=args.jobs,
            version=args.gdc_client_version,
            install_dir=args.install_dir,
        )
        print_reports(cohort)
    elif args.command == "cases":
        cohort.cases()
    elif args.command == "clinical":
        cohort.clinical()
    elif args.command == "files":
        cohort.maf_files()
        cohort.duplicates()
    elif args.command == "select":
        cohort.selection()
    elif args.command == "download":
        cohort.download(
            gdc_client=args.gdc_client,
            n_clients=args.jobs,
            download=not args.no_download,
            version=args.gdc_client_version,
            install_dir=args.install_dir,
        )
        cohort.write_download_record(record_dir=args.record_dir)
    elif args.command == "merge":
        cohort.merged_maf()
        cohort.samples_for_review()
        cohort.write_download_record(record_dir=getattr(args, "record_dir", None))
    elif args.command == "reports":
        print_reports(cohort)

    if not args.no_log:
        logger.info(f"\nLog: {cohort.log_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
