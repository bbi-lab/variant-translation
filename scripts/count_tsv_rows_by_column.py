"""Count rows in a TSV file, grouped by the value in a specified column.

Output is a two-column TSV with the group value and row count.  By default
results are sorted alphabetically by value and written to stdout.

Example:
  python scripts/count_tsv_rows_by_column.py \\
    --input reverse_translated_variants.tsv \\
    --column Gene \\
    --output gene_counts.tsv
"""

from __future__ import annotations

import csv
import sys
from collections import Counter
from contextlib import contextmanager
from pathlib import Path
from typing import Generator, TextIO

import click


@contextmanager
def open_output(output_path: Path | None) -> Generator[TextIO, None, None]:
    if output_path is None:
        yield sys.stdout
    else:
        with output_path.open("w", newline="") as handle:
            yield handle


SORT_CHOICES = click.Choice(
    ["value-asc", "value-desc", "count-asc", "count-desc"],
    case_sensitive=False,
)


def configure_csv_field_size_limit(field_size_limit: int) -> None:
    try:
        csv.field_size_limit(field_size_limit)
    except OverflowError as exception:
        raise click.ClickException(
            f"Invalid --csv-field-size-limit value {field_size_limit}: exceeds platform maximum."
        ) from exception


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
    help="Name of the column to group by.",
)
@click.option(
    "--output",
    "output_path",
    default=None,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Output TSV file path.  Defaults to stdout.",
)
@click.option(
    "--csv-field-size-limit",
    default=csv.field_size_limit(),
    show_default=True,
    type=click.IntRange(min=1),
    help="Maximum per-field character length for CSV/TSV parsing in batch mode.",
)
@click.option(
    "--sort",
    "sort_order",
    default="value-asc",
    show_default=True,
    type=SORT_CHOICES,
    help=(
        "Sort order for the output rows.  "
        "value-asc / value-desc sort alphabetically by the group value; "
        "count-asc / count-desc sort numerically by row count."
    ),
)
def main(
    input_path: Path,
    column_name: str,
    csv_field_size_limit: int,
    output_path: Path | None,
    sort_order: str,
) -> None:
    configure_csv_field_size_limit(csv_field_size_limit)

    with input_path.open("r", newline="") as input_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        if reader.fieldnames is None:
            raise click.ClickException("Input TSV has no header row.")

        if column_name not in reader.fieldnames:
            raise click.ClickException(f"Column not found in input TSV: {column_name}")

        counts: Counter[str] = Counter()
        total_rows = 0
        for row in reader:
            total_rows += 1
            counts[row[column_name]] += 1

    # Sort
    sort_key, sort_reverse = sort_order.split("-")
    reverse = sort_reverse == "desc"

    if sort_key == "value":
        sorted_items = sorted(counts.items(), key=lambda item: item[0], reverse=reverse)
    else:  # count
        sorted_items = sorted(counts.items(), key=lambda item: item[1], reverse=reverse)

    with open_output(output_path) as out_handle:
        writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
        writer.writerow([column_name, "count"])
        for value, count in sorted_items:
            writer.writerow([value, count])

    distinct_values = len(counts)
    click.echo(
        f"{total_rows} input rows, {distinct_values} distinct values in column '{column_name}'.",
        err=True,
    )


if __name__ == "__main__":
    main()
