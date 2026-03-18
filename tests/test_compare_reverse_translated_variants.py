import csv
from pathlib import Path

from click.testing import CliRunner
import pytest

from src.scripts import compare_reverse_translated_variants as cmp


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


def write_delimited(path: Path, rows: list[dict[str, str]], fieldnames: list[str], delimiter: str) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        assert reader.fieldnames is not None
        return list(reader.fieldnames), rows


def test_compare_outputs_and_pass_through_with_reference_ignoring(runner: CliRunner, tmp_path: Path) -> None:
    a_input = tmp_path / "a.tsv"
    b_input = tmp_path / "b.tsv"
    a_only_output = tmp_path / "a_only.tsv"
    b_only_output = tmp_path / "b_only.tsv"
    different_output = tmp_path / "different.tsv"

    write_delimited(
        a_input,
        rows=[
            {
                "sample_id": "a1",
                "note": "equivalent",
                "transcript": "NM_000001.1",
                "hgvs_p": "p.Arg1His",
                "hgvs_c": "NM_000001.1:c.1A>G|NM_000001.1:c.1A>T",
                "hgvs_g": "NC_000001.11:g.1A>G|NC_000001.11:g.1A>T",
            },
            {
                "sample_id": "a2",
                "note": "only a",
                "transcript": "NM_000002.1",
                "hgvs_p": "p.Arg2His",
                "hgvs_c": "NM_000002.1:c.2A>G",
                "hgvs_g": "NC_000001.11:g.2A>G",
            },
            {
                "sample_id": "a3",
                "note": "different",
                "transcript": "NM_000003.1",
                "hgvs_p": "p.Arg3His",
                "hgvs_c": "NM_000003.1:c.3A>G",
                "hgvs_g": "NC_000001.11:g.3A>G",
            },
        ],
        fieldnames=["sample_id", "note", "transcript", "hgvs_p", "hgvs_c", "hgvs_g"],
        delimiter="\t",
    )

    write_delimited(
        b_input,
        rows=[
            {
                "record_id": "b1",
                "group": "equivalent",
                "transcript": "NM_000001.1",
                "hgvs_p": "p.Arg1His",
                "hgvs_c": "c.1A>G|c.1A>T",
                "hgvs_g": "g.1A>G|g.1A>T",
            },
            {
                "record_id": "b2",
                "group": "only b",
                "transcript": "NM_000004.1",
                "hgvs_p": "p.Arg4His",
                "hgvs_c": "c.4A>G",
                "hgvs_g": "g.4A>G",
            },
            {
                "record_id": "b3",
                "group": "different",
                "transcript": "NM_000003.1",
                "hgvs_p": "p.Arg3His",
                "hgvs_c": "c.3A>C",
                "hgvs_g": "g.3A>G",
            },
        ],
        fieldnames=["record_id", "group", "transcript", "hgvs_p", "hgvs_c", "hgvs_g"],
        delimiter="\t",
    )

    result = runner.invoke(
        cmp.main,
        [
            "--a-input",
            str(a_input),
            "--b-input",
            str(b_input),
            "--ignore-reference-differences",
            "--pass-through-all-columns",
            "--a-prefix",
            "a_",
            "--b-prefix",
            "b_",
            "--a-pass-through-prefix",
            "a_extra_",
            "--b-pass-through-prefix",
            "b_extra_",
            "--a-only-output",
            str(a_only_output),
            "--b-only-output",
            str(b_only_output),
            "--different-output",
            str(different_output),
        ],
    )

    assert result.exit_code == 0, result.output

    a_only_fields, a_only_rows = read_tsv(a_only_output)
    b_only_fields, b_only_rows = read_tsv(b_only_output)
    diff_fields, diff_rows = read_tsv(different_output)

    assert "sample_id" in a_only_fields
    assert len(a_only_rows) == 1
    assert a_only_rows[0]["transcript"] == "NM_000002.1"

    assert "record_id" in b_only_fields
    assert len(b_only_rows) == 1
    assert b_only_rows[0]["transcript"] == "NM_000004.1"

    assert len(diff_rows) == 1
    assert "a_transcript" in diff_fields
    assert "a_hgvs_p" in diff_fields
    assert "a_hgvs_c" in diff_fields
    assert "a_hgvs_g" in diff_fields
    assert "b_transcript" in diff_fields
    assert "b_hgvs_p" in diff_fields
    assert "b_hgvs_c" in diff_fields
    assert "b_hgvs_g" in diff_fields
    assert "a_extra_sample_id" in diff_fields
    assert "b_extra_record_id" in diff_fields

    assert diff_rows[0]["a_transcript"] == "NM_000003.1"
    assert diff_rows[0]["b_transcript"] == "NM_000003.1"
    assert diff_rows[0]["a_hgvs_c"] == "NM_000003.1:c.3A>G"
    assert diff_rows[0]["b_hgvs_c"] == "c.3A>C"


def test_compare_hgvs_c_only_does_not_require_hgvs_g_columns(runner: CliRunner, tmp_path: Path) -> None:
    a_input = tmp_path / "a.csv"
    b_input = tmp_path / "b.csv"
    a_only_output = tmp_path / "a_only.tsv"
    b_only_output = tmp_path / "b_only.tsv"
    different_output = tmp_path / "different.tsv"

    write_delimited(
        a_input,
        rows=[
            {
                "transcript": "NM_000010.1",
                "hgvs_p": "p.Arg10His",
                "hgvs_c": "c.10A>G",
            }
        ],
        fieldnames=["transcript", "hgvs_p", "hgvs_c"],
        delimiter=",",
    )

    write_delimited(
        b_input,
        rows=[
            {
                "transcript": "NM_000010.1",
                "hgvs_p": "p.Arg10His",
                "hgvs_c": "c.10A>T",
            }
        ],
        fieldnames=["transcript", "hgvs_p", "hgvs_c"],
        delimiter=",",
    )

    result = runner.invoke(
        cmp.main,
        [
            "--a-input",
            str(a_input),
            "--a-input-format",
            "csv",
            "--b-input",
            str(b_input),
            "--b-input-format",
            "csv",
            "--compare",
            "hgvs_c",
            "--a-only-output",
            str(a_only_output),
            "--b-only-output",
            str(b_only_output),
            "--different-output",
            str(different_output),
        ],
    )

    assert result.exit_code == 0, result.output

    _, a_only_rows = read_tsv(a_only_output)
    _, b_only_rows = read_tsv(b_only_output)
    diff_fields, diff_rows = read_tsv(different_output)

    assert len(a_only_rows) == 0
    assert len(b_only_rows) == 0
    assert len(diff_rows) == 1
    assert "a_hgvs_g" not in diff_fields
    assert "b_hgvs_g" not in diff_fields
    assert diff_rows[0]["a_hgvs_c"] == "c.10A>G"
    assert diff_rows[0]["b_hgvs_c"] == "c.10A>T"


def test_compare_pass_through_columns_accepts_csv_list_and_repeats(runner: CliRunner, tmp_path: Path) -> None:
    a_input = tmp_path / "a.tsv"
    b_input = tmp_path / "b.tsv"
    a_only_output = tmp_path / "a_only.tsv"
    b_only_output = tmp_path / "b_only.tsv"
    different_output = tmp_path / "different.tsv"

    write_delimited(
        a_input,
        rows=[
            {
                "dataset": "set1",
                "sample_id": "a1",
                "batch": "ba",
                "transcript": "NM_000020.1",
                "hgvs_p": "p.Arg20His",
                "hgvs_c": "c.20A>G",
                "hgvs_g": "g.20A>G",
            }
        ],
        fieldnames=["dataset", "sample_id", "batch", "transcript", "hgvs_p", "hgvs_c", "hgvs_g"],
        delimiter="\t",
    )

    write_delimited(
        b_input,
        rows=[
            {
                "dataset": "set1",
                "sample_id": "b1",
                "batch": "bb",
                "transcript": "NM_000020.1",
                "hgvs_p": "p.Arg20His",
                "hgvs_c": "c.20A>T",
                "hgvs_g": "g.20A>G",
            }
        ],
        fieldnames=["dataset", "sample_id", "batch", "transcript", "hgvs_p", "hgvs_c", "hgvs_g"],
        delimiter="\t",
    )

    result = runner.invoke(
        cmp.main,
        [
            "--a-input",
            str(a_input),
            "--b-input",
            str(b_input),
            "--compare",
            "hgvs_c",
            "--pass-through-columns",
            "dataset,sample_id",
            "--pass-through-columns",
            "batch",
            "--a-pass-through-prefix",
            "a_extra_",
            "--b-pass-through-prefix",
            "b_extra_",
            "--a-only-output",
            str(a_only_output),
            "--b-only-output",
            str(b_only_output),
            "--different-output",
            str(different_output),
        ],
    )

    assert result.exit_code == 0, result.output

    diff_fields, diff_rows = read_tsv(different_output)
    assert len(diff_rows) == 1
    assert "a_extra_dataset" in diff_fields
    assert "a_extra_sample_id" in diff_fields
    assert "a_extra_batch" in diff_fields
    assert "b_extra_dataset" in diff_fields
    assert "b_extra_sample_id" in diff_fields
    assert "b_extra_batch" in diff_fields
    assert diff_rows[0]["a_extra_dataset"] == "set1"
    assert diff_rows[0]["b_extra_dataset"] == "set1"


def test_compare_match_columns_disambiguates_overlapping_variants(runner: CliRunner, tmp_path: Path) -> None:
    a_input = tmp_path / "a.tsv"
    b_input = tmp_path / "b.tsv"
    a_only_output = tmp_path / "a_only.tsv"
    b_only_output = tmp_path / "b_only.tsv"
    different_output = tmp_path / "different.tsv"

    write_delimited(
        a_input,
        rows=[
            {"dataset": "set1", "transcript": "NM_000030.1", "hgvs_p": "p.Arg30His", "hgvs_c": "c.30A>G"},
            {"dataset": "set2", "transcript": "NM_000030.1", "hgvs_p": "p.Arg30His", "hgvs_c": "c.30A>T"},
        ],
        fieldnames=["dataset", "transcript", "hgvs_p", "hgvs_c"],
        delimiter="\t",
    )

    write_delimited(
        b_input,
        rows=[
            {"dataset": "set1", "transcript": "NM_000030.1", "hgvs_p": "p.Arg30His", "hgvs_c": "c.30A>G"},
            {"dataset": "set2", "transcript": "NM_000030.1", "hgvs_p": "p.Arg30His", "hgvs_c": "c.30A>C"},
        ],
        fieldnames=["dataset", "transcript", "hgvs_p", "hgvs_c"],
        delimiter="\t",
    )

    result = runner.invoke(
        cmp.main,
        [
            "--a-input",
            str(a_input),
            "--b-input",
            str(b_input),
            "--compare",
            "hgvs_c",
            "--match-columns",
            "dataset",
            "--pass-through-columns",
            "dataset",
            "--a-only-output",
            str(a_only_output),
            "--b-only-output",
            str(b_only_output),
            "--different-output",
            str(different_output),
        ],
    )

    assert result.exit_code == 0, result.output

    _, a_only_rows = read_tsv(a_only_output)
    _, b_only_rows = read_tsv(b_only_output)
    diff_fields, diff_rows = read_tsv(different_output)

    assert len(a_only_rows) == 0
    assert len(b_only_rows) == 0
    assert len(diff_rows) == 1
    assert "a_extra_dataset" in diff_fields
    assert "b_extra_dataset" in diff_fields
    assert diff_rows[0]["a_extra_dataset"] == "set2"
    assert diff_rows[0]["a_hgvs_c"] == "c.30A>T"
    assert diff_rows[0]["b_hgvs_c"] == "c.30A>C"


def test_compare_handles_large_csv_field_values(runner: CliRunner, tmp_path: Path) -> None:
    a_input = tmp_path / "a.tsv"
    b_input = tmp_path / "b.tsv"
    a_only_output = tmp_path / "a_only.tsv"
    b_only_output = tmp_path / "b_only.tsv"
    different_output = tmp_path / "different.tsv"

    very_large_variant = "c." + ("A" * 200000)

    write_delimited(
        a_input,
        rows=[
            {
                "transcript": "NM_999999.1",
                "hgvs_p": "p.Arg1His",
                "hgvs_c": very_large_variant,
            }
        ],
        fieldnames=["transcript", "hgvs_p", "hgvs_c"],
        delimiter="\t",
    )

    write_delimited(
        b_input,
        rows=[
            {
                "transcript": "NM_999999.1",
                "hgvs_p": "p.Arg1His",
                "hgvs_c": very_large_variant,
            }
        ],
        fieldnames=["transcript", "hgvs_p", "hgvs_c"],
        delimiter="\t",
    )

    result = runner.invoke(
        cmp.main,
        [
            "--a-input",
            str(a_input),
            "--b-input",
            str(b_input),
            "--compare",
            "hgvs_c",
            "--a-only-output",
            str(a_only_output),
            "--b-only-output",
            str(b_only_output),
            "--different-output",
            str(different_output),
        ],
    )

    assert result.exit_code == 0, result.output

    _, a_only_rows = read_tsv(a_only_output)
    _, b_only_rows = read_tsv(b_only_output)
    _, diff_rows = read_tsv(different_output)

    assert len(a_only_rows) == 0
    assert len(b_only_rows) == 0
    assert len(diff_rows) == 0


def test_compare_writes_missing_transcript_or_hgvs_p_rows(runner: CliRunner, tmp_path: Path) -> None:
    a_input = tmp_path / "a.tsv"
    b_input = tmp_path / "b.tsv"
    a_only_output = tmp_path / "a_only.tsv"
    b_only_output = tmp_path / "b_only.tsv"
    different_output = tmp_path / "different.tsv"
    a_missing_output = tmp_path / "a_missing.tsv"
    b_missing_output = tmp_path / "b_missing.tsv"

    write_delimited(
        a_input,
        rows=[
            {"transcript": "NM_000040.1", "hgvs_p": "p.Arg40His", "hgvs_c": "c.40A>G"},
            {"transcript": "", "hgvs_p": "p.Arg41His", "hgvs_c": "c.41A>G"},
            {"transcript": "NM_000042.1", "hgvs_p": "", "hgvs_c": "c.42A>G"},
        ],
        fieldnames=["transcript", "hgvs_p", "hgvs_c"],
        delimiter="\t",
    )

    write_delimited(
        b_input,
        rows=[
            {"transcript": "NM_000040.1", "hgvs_p": "p.Arg40His", "hgvs_c": "c.40A>G"},
            {"transcript": "", "hgvs_p": "p.Arg43His", "hgvs_c": "c.43A>G"},
        ],
        fieldnames=["transcript", "hgvs_p", "hgvs_c"],
        delimiter="\t",
    )

    result = runner.invoke(
        cmp.main,
        [
            "--a-input",
            str(a_input),
            "--b-input",
            str(b_input),
            "--compare",
            "hgvs_c",
            "--a-only-output",
            str(a_only_output),
            "--b-only-output",
            str(b_only_output),
            "--different-output",
            str(different_output),
            "--a-missing-output",
            str(a_missing_output),
            "--b-missing-output",
            str(b_missing_output),
        ],
    )

    assert result.exit_code == 0, result.output

    _, a_only_rows = read_tsv(a_only_output)
    _, b_only_rows = read_tsv(b_only_output)
    _, diff_rows = read_tsv(different_output)
    _, a_missing_rows = read_tsv(a_missing_output)
    _, b_missing_rows = read_tsv(b_missing_output)

    assert len(a_only_rows) == 0
    assert len(b_only_rows) == 0
    assert len(diff_rows) == 0
    assert len(a_missing_rows) == 2
    assert len(b_missing_rows) == 1
