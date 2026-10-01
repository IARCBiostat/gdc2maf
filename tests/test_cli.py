"""The gdc2maf command line: parsing, defaults, and the cohort it builds."""

import contextlib
import io

import pytest

from gdc2maf.cli import build_parser, cohort_from_args, default_cohort_name


def parse(*argv):
    """Parse a command line."""
    return build_parser().parse_args(list(argv))


def test_every_documented_subcommand_exists():
    commands = ["maf", "cases", "clinical", "files", "select", "download",
                "merge", "reports", "client"]
    for command in commands:
        extra = [] if command == "client" else ["--project", "TCGA-TEST"]
        assert parse(command, *extra)


def test_a_command_is_required():
    with pytest.raises(SystemExit):
        parse()


def test_a_project_is_required_for_cohort_commands():
    with pytest.raises(SystemExit):
        parse("cases")


def test_projects_accumulate():
    args = parse("cases", "--project", "TCGA-LUAD", "--project", "TCGA-LUSC")
    assert args.project == ["TCGA-LUAD", "TCGA-LUSC"]


def test_an_invalid_sex_is_rejected():
    with pytest.raises(SystemExit):
        parse("cases", "--project", "TCGA-TEST", "--sex-at-birth", "other")


def test_an_invalid_refresh_step_is_rejected_by_the_parser():
    with pytest.raises(SystemExit):
        parse("cases", "--project", "TCGA-TEST", "--refresh", "nonsense")


def test_default_cohort_name_joins_projects_and_sex():
    assert default_cohort_name(["TCGA-LUAD", "TCGA-LUSC"], ["male"]) == (
        "TCGA-LUAD_TCGA-LUSC_MALE"
    )
    assert default_cohort_name(["TCGA-BRCA"], None) == "TCGA-BRCA"


def test_the_cohort_output_dir_is_the_parent_plus_the_name(tmp_path):
    args = parse(
        "maf", "--project", "TCGA-LUAD", "--sex-at-birth", "male",
        "--name", "Lung_MALE", "--out", str(tmp_path),
    )
    cohort = cohort_from_args(args)
    assert cohort.name == "Lung_MALE"
    assert cohort.out_dir == f"{tmp_path}/Lung_MALE"


def test_cohort_options_reach_the_cohort(tmp_path):
    args = parse(
        "maf", "--project", "TCGA-LUAD", "--out", str(tmp_path),
        "--sex-fallback", "pancan.tsv",
        "--remove-gdc-filter-flag", "NonExonic",
        "--remove-gdc-filter-flag", "gdc_pon",
        "--refresh", "files",
    )
    cohort = cohort_from_args(args)
    assert cohort.sex_fallback_path == "pancan.tsv"
    assert cohort.remove_gdc_filter_flags == ["NonExonic", "gdc_pon"]
    assert cohort.refresh == ["files"]


def test_download_defaults(tmp_path):
    args = parse("download", "--project", "TCGA-TEST", "--out", str(tmp_path))
    assert args.jobs == 8
    assert args.gdc_client_version == "2.3"
    assert args.no_download is False


def test_omitting_sex_keeps_every_sex(tmp_path):
    args = parse("maf", "--project", "TCGA-BRCA", "--out", str(tmp_path))
    assert args.sex_at_birth is None
    assert cohort_from_args(args).sex_at_birth is None


def test_the_help_states_the_all_sexes_default():
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer), pytest.raises(SystemExit):
        build_parser().parse_args(["cases", "--help"])
    # argparse wraps to the terminal width, so compare on collapsed whitespace
    printed = " ".join(buffer.getvalue().split())
    assert "every sex is kept" in printed


def test_help_is_available_as_both_h_and_help():
    for flag in ["-h", "--help"]:
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer), pytest.raises(SystemExit) as exit_info:
            build_parser().parse_args([flag])
        assert exit_info.value.code == 0
        assert "usage: gdc2maf" in buffer.getvalue()


def test_help_is_available_for_every_subcommand():
    commands = ["maf", "cases", "clinical", "files", "select", "download",
                "merge", "reports", "client"]
    for command in commands:
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer), pytest.raises(SystemExit) as exit_info:
            build_parser().parse_args([command, "-h"])
        assert exit_info.value.code == 0
        assert f"usage: gdc2maf {command}" in buffer.getvalue()


def help_text(*argv):
    """Capture the help printed for a command line."""
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer), pytest.raises(SystemExit):
        build_parser().parse_args(list(argv))
    return buffer.getvalue()


def test_help_explains_where_to_find_project_ids():
    printed = help_text("cases", "-h")
    assert "finding project IDs" in printed
    assert "https://portal.gdc.cancer.gov/projects" in printed
    assert "https://api.gdc.cancer.gov/projects" in printed


def test_help_warns_that_a_primary_site_is_not_a_disease():
    printed = " ".join(help_text("maf", "-h").split())
    assert "a primary site is not a disease" in printed
    assert "TCGA-LUAD" in printed and "TCGA-LUSC" in printed


def test_the_example_api_query_survives_copy_and_paste():
    # the filters argument must stay on one line, and size/format must be
    # separate --data-urlencode values (combining them with & encodes the &)
    printed = help_text("cases", "-h")
    filter_lines = [ln for ln in printed.splitlines() if "primary_site" in ln]
    assert len(filter_lines) == 1
    assert filter_lines[0].rstrip().endswith("}}'")
    assert "size=100&format" not in printed


def test_project_guidance_is_on_every_command_that_takes_projects():
    for command in ["maf", "cases", "clinical", "files", "select", "download",
                    "merge", "reports"]:
        assert "finding project IDs" in help_text(command, "-h"), command


def test_the_client_command_does_not_carry_project_guidance():
    # it takes no --project, so the guidance would be noise
    assert "finding project IDs" not in help_text("client", "-h")
