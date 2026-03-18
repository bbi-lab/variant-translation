"""Split a TSV file into three files by variant count in a selected column.

Rows are classified as:
- empty: the selected column is missing/blank after trimming
- single: exactly one non-empty variant token
- multiple: two or more non-empty variant tokens separated by the given separator

Example:
  python scripts/split_tsv_by_variant_count.py \
    --input reverse_translated_variants.tsv \
    --column hgvs_c \
    --separator "^" \
    --empty-output hgvs_c_empty.tsv \
    --single-output hgvs_c_single.tsv \
    --multiple-output hgvs_c_multiple.tsv
"""

from __future__ import annotations

import csv
from pathlib import Path

import click


def classify_variant_cell(cell_value: str | None, separator: str) -> str:
    normalized_value = (cell_value or "").strip()
    if not normalized_value:
        return "empty"

    tokens = [token.strip() for token in normalized_value.split(separator) if token.strip()]
    if len(tokens) <= 1:
        return "single"

    return "multiple"


@click.command()
@click.option(
    "--input",
    "input_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Input TSV file path.",
)
@click.option(
    "--column",
    "column_name",
    required=True,
    type=str,
    help="Name of the TSV column that contains joined variant strings.",
)
@click.option(
    "--separator",
    default="^",
    show_default=True,
    type=str,
    help="Single-character separator between variant strings in the selected column.",
)
@click.option(
    "--empty-output",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Output TSV for rows where the selected column is blank.",
)
@click.option(
    "--single-output",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Output TSV for rows where the selected column has one variant string.",
)
@click.option(
    "--multiple-output",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Output TSV for rows where the selected column has multiple variant strings.",
)
def main(
    input_path: Path,
    column_name: str,
    separator: str,
    empty_output: Path,
    single_output: Path,
    multiple_output: Path,
) -> None:
    if len(separator) != 1:
        raise click.ClickException("--separator must be exactly one character.")

    output_paths = [empty_output, single_output, multiple_output]
    if len({str(path) for path in output_paths}) != 3:
        raise click.ClickException("--empty-output, --single-output, and --multiple-output must be different paths.")

    with input_path.open("r", newline="") as input_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        if reader.fieldnames is None:
            raise click.ClickException("Input TSV has no header row.")

        field_names = list(reader.fieldnames)
        if column_name not in field_names:
            raise click.ClickException(f"Column not found in input TSV: {column_name}")

        with (
            empty_output.open("w", newline="") as empty_handle,
            single_output.open("w", newline="") as single_handle,
            multiple_output.open("w", newline="") as multiple_handle,
        ):
            empty_writer = csv.DictWriter(empty_handle, fieldnames=field_names, delimiter="\t")
            single_writer = csv.DictWriter(single_handle, fieldnames=field_names, delimiter="\t")
            multiple_writer = csv.DictWriter(multiple_handle, fieldnames=field_names, delimiter="\t")

            empty_writer.writeheader()
            single_writer.writeheader()
            multiple_writer.writeheader()

            total_rows = 0
            empty_rows = 0
            single_rows = 0
            multiple_rows = 0

            for row in reader:
                total_rows += 1
                classification = classify_variant_cell(row.get(column_name), separator=separator)

                if classification == "empty":
                    empty_writer.writerow(row)
                    empty_rows += 1
                elif classification == "single":
                    single_writer.writerow(row)
                    single_rows += 1
                else:
                    multiple_writer.writerow(row)
                    multiple_rows += 1

    click.echo(
        "Split summary: "
        f"{total_rows} input rows, "
        f"{empty_rows} empty, "
        f"{single_rows} single, "
        f"{multiple_rows} multiple.",
        err=True,
    )


if __name__ == "__main__":
    main()
