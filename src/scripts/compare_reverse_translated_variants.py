"""
Compare two files containing reverse-translated protein variants.

This script compares single-line reverse-translation outputs from two files (A and B),
matching rows by exact transcript + HGVS p. values (and optional additional key columns) and reporting:
- rows present only in A,
- rows present only in B,
- rows present in both but with different reverse translations.

Examples:
  python -m src.scripts.compare_reverse_translated_variants \
    --a-input a.tsv \
    --b-input b.tsv \
    --a-only-output a_only.tsv \
    --b-only-output b_only.tsv \
    --different-output differences.tsv

  python -m src.scripts.compare_reverse_translated_variants \
    --a-input a.csv --a-input-format csv --a-separator "^" \
    --b-input b.tsv --b-input-format tsv --b-separator "|" \
    --compare both \
    --ignore-reference-differences \
	--pass-through-all-columns \
    --a-pass-through-prefix a_extra_ \
    --b-pass-through-prefix b_extra_ \
    --a-only-output a_only.tsv \
    --b-only-output b_only.tsv \
    --different-output differences.tsv
"""

import csv
from dataclasses import dataclass
from pathlib import Path
import sys

import click


@dataclass(frozen=True)
class InputSpec:
	label: str
	path: Path
	input_format: str
	separator: str
	transcript_column: str
	hgvs_p_column: str
	hgvs_c_column: str
	hgvs_g_column: str


def normalize_hgvs_reference(token: str, molecule_type: str, ignore_reference_differences: bool) -> str:
	normalized_token = token.strip()
	if not ignore_reference_differences:
		return normalized_token

	if ":" not in normalized_token:
		return normalized_token

	_, suffix = normalized_token.split(":", 1)
	if suffix.startswith(f"{molecule_type}."):
		return suffix

	return normalized_token


def split_joined_variants(
	joined_variants: str,
	separator: str,
	molecule_type: str,
	ignore_reference_differences: bool,
) -> set[str]:
	normalized_joined_variants = (joined_variants or "").strip()
	if not normalized_joined_variants:
		return set()

	return {
		normalize_hgvs_reference(token, molecule_type=molecule_type, ignore_reference_differences=ignore_reference_differences)
		for token in normalized_joined_variants.split(separator)
		if token.strip()
	}


def ensure_required_columns(field_names: list[str], required_columns: list[str], label: str) -> None:
	missing_columns = [column for column in required_columns if column not in field_names]
	if missing_columns:
		raise click.ClickException(
			f"Missing required columns in {label}: {', '.join(missing_columns)}"
		)


def parse_column_list_option(raw_values: tuple[str, ...], option_name: str) -> list[str]:
	parsed_values: list[str] = []
	seen_values: set[str] = set()

	for raw_group in raw_values:
		for raw_value in raw_group.split(","):
			value = raw_value.strip()
			if not value:
				raise click.ClickException(
					f"{option_name} entries cannot be empty. Use a comma-separated list like 'col1,col2'."
				)
			if value in seen_values:
				continue
			seen_values.add(value)
			parsed_values.append(value)

	return parsed_values


def configure_csv_field_size_limit() -> None:
	field_size_limit = sys.maxsize
	while True:
		try:
			csv.field_size_limit(field_size_limit)
			return
		except OverflowError:
			field_size_limit //= 10
			if field_size_limit < 1:
				raise click.ClickException("Unable to configure CSV field size limit for large fields.")


def read_rows_by_key(
	spec: InputSpec,
	compare: str,
	match_columns: list[str],
) -> tuple[list[str], dict[tuple[str, ...], dict[str, str]], list[tuple[str, ...]], list[dict[str, str]]]:
	input_delimiter = "\t" if spec.input_format == "tsv" else ","

	with spec.path.open("r", newline="") as input_handle:
		reader = csv.DictReader(input_handle, delimiter=input_delimiter)
		if reader.fieldnames is None:
			raise click.ClickException(f"{spec.label} input has no header row.")

		field_names = list(reader.fieldnames)
		required_columns = [spec.transcript_column, spec.hgvs_p_column]
		required_columns.extend(match_columns)
		if compare in {"hgvs_c", "both"}:
			required_columns.append(spec.hgvs_c_column)
		if compare in {"hgvs_g", "both"}:
			required_columns.append(spec.hgvs_g_column)
		ensure_required_columns(field_names, required_columns, label=spec.label)

		rows_by_key: dict[tuple[str, ...], dict[str, str]] = {}
		row_key_order: list[tuple[str, ...]] = []
		missing_key_rows: list[dict[str, str]] = []
		key_column_names = [spec.transcript_column, spec.hgvs_p_column, *match_columns]

		for row_number, row in enumerate(reader, start=2):
			transcript_value = (row.get(spec.transcript_column) or "").strip()
			hgvs_p_value = (row.get(spec.hgvs_p_column) or "").strip()

			if not transcript_value or not hgvs_p_value:
				missing_key_rows.append(dict(row))
				continue

			row_key_values = [transcript_value, hgvs_p_value]
			row_key_values.extend((row.get(column_name) or "").strip() for column_name in match_columns)
			row_key = tuple(row_key_values)
			if row_key in rows_by_key:
				formatted_key = ", ".join(
					f"{column_name}={value}" for column_name, value in zip(key_column_names, row_key_values)
				)
				raise click.ClickException(
					f"Duplicate match key in {spec.label} at row {row_number}: {formatted_key}"
				)

			rows_by_key[row_key] = row
			row_key_order.append(row_key)

	return field_names, rows_by_key, row_key_order, missing_key_rows


def rows_differ(
	a_row: dict[str, str],
	b_row: dict[str, str],
	a_spec: InputSpec,
	b_spec: InputSpec,
	compare: str,
	ignore_reference_differences: bool,
) -> bool:
	if compare in {"hgvs_c", "both"}:
		a_hgvs_c = split_joined_variants(
			a_row.get(a_spec.hgvs_c_column, ""),
			separator=a_spec.separator,
			molecule_type="c",
			ignore_reference_differences=ignore_reference_differences,
		)
		b_hgvs_c = split_joined_variants(
			b_row.get(b_spec.hgvs_c_column, ""),
			separator=b_spec.separator,
			molecule_type="c",
			ignore_reference_differences=ignore_reference_differences,
		)
		if a_hgvs_c != b_hgvs_c:
			return True

	if compare in {"hgvs_g", "both"}:
		a_hgvs_g = split_joined_variants(
			a_row.get(a_spec.hgvs_g_column, ""),
			separator=a_spec.separator,
			molecule_type="g",
			ignore_reference_differences=ignore_reference_differences,
		)
		b_hgvs_g = split_joined_variants(
			b_row.get(b_spec.hgvs_g_column, ""),
			separator=b_spec.separator,
			molecule_type="g",
			ignore_reference_differences=ignore_reference_differences,
		)
		if a_hgvs_g != b_hgvs_g:
			return True

	return False


def build_prefixed_column_names(columns: list[str], prefix: str) -> dict[str, str]:
	return {column: f"{prefix}{column}" for column in columns}


def write_tsv(path: Path, field_names: list[str], rows: list[dict[str, str]]) -> None:
	with path.open("w", newline="") as output_handle:
		writer = csv.DictWriter(output_handle, fieldnames=field_names, delimiter="\t")
		writer.writeheader()
		writer.writerows(rows)


@click.command()
@click.option("--a-input", "a_input", required=True, type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("--b-input", "b_input", required=True, type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option(
	"--a-input-format",
	default="tsv",
	show_default=True,
	type=click.Choice(["tsv", "csv"], case_sensitive=False),
)
@click.option(
	"--b-input-format",
	default="tsv",
	show_default=True,
	type=click.Choice(["tsv", "csv"], case_sensitive=False),
)
@click.option("--a-separator", default="|", show_default=True, type=str, help="Joined reverse-translation separator in A.")
@click.option("--b-separator", default="|", show_default=True, type=str, help="Joined reverse-translation separator in B.")
@click.option("--a-transcript-column", default="transcript", show_default=True, type=str)
@click.option("--b-transcript-column", default="transcript", show_default=True, type=str)
@click.option("--a-hgvs-p-column", default="hgvs_p", show_default=True, type=str)
@click.option("--b-hgvs-p-column", default="hgvs_p", show_default=True, type=str)
@click.option("--a-hgvs-c-column", default="hgvs_c", show_default=True, type=str)
@click.option("--b-hgvs-c-column", default="hgvs_c", show_default=True, type=str)
@click.option("--a-hgvs-g-column", default="hgvs_g", show_default=True, type=str)
@click.option("--b-hgvs-g-column", default="hgvs_g", show_default=True, type=str)
@click.option(
	"--compare",
	default="both",
	show_default=True,
	type=click.Choice(["hgvs_c", "hgvs_g", "both"], case_sensitive=False),
	help="Which reverse-translation columns to compare.",
)
@click.option(
	"--match-columns",
	"match_columns",
	multiple=True,
	type=str,
	help=(
		"Additional columns that must match between A and B. "
		"Accepts comma-separated values and may be repeated."
	),
)
@click.option(
	"--ignore-reference-differences",
	is_flag=True,
	default=False,
	help="Ignore accession/reference prefixes (e.g., NM_:c. vs c., NC_:g. vs g.) when comparing variants.",
)
@click.option("--a-prefix", default="a_", show_default=True, type=str, help="Prefix for A core columns in differences output.")
@click.option("--b-prefix", default="b_", show_default=True, type=str, help="Prefix for B core columns in differences output.")
@click.option(
	"--pass-through-columns",
	"pass_through_selected_columns",
	multiple=True,
	type=str,
	help=(
		"Include selected non-core columns from A and/or B in differences output. "
		"Accepts comma-separated values and may be repeated."
	),
)
@click.option(
	"--pass-through-all-columns",
	is_flag=True,
	default=False,
	help="Include all non-core columns from A and B in differences output.",
)
@click.option(
	"--a-pass-through-prefix",
	default="a_extra_",
	show_default=True,
	type=str,
	help="Prefix for A non-core columns in differences output when --pass-through-columns is enabled.",
)
@click.option(
	"--b-pass-through-prefix",
	default="b_extra_",
	show_default=True,
	type=str,
	help="Prefix for B non-core columns in differences output when --pass-through-columns is enabled.",
)
@click.option("--a-only-output", required=True, type=click.Path(dir_okay=False, path_type=Path))
@click.option("--b-only-output", required=True, type=click.Path(dir_okay=False, path_type=Path))
@click.option("--different-output", required=True, type=click.Path(dir_okay=False, path_type=Path))
@click.option(
	"--a-missing-output",
	default="a_missing_transcript_or_hgvs_p.tsv",
	show_default=True,
	type=click.Path(dir_okay=False, path_type=Path),
	help="Output TSV for A rows missing transcript or hgvs_p.",
)
@click.option(
	"--b-missing-output",
	default="b_missing_transcript_or_hgvs_p.tsv",
	show_default=True,
	type=click.Path(dir_okay=False, path_type=Path),
	help="Output TSV for B rows missing transcript or hgvs_p.",
)
def main(
	a_input: Path,
	b_input: Path,
	a_input_format: str,
	b_input_format: str,
	a_separator: str,
	b_separator: str,
	a_transcript_column: str,
	b_transcript_column: str,
	a_hgvs_p_column: str,
	b_hgvs_p_column: str,
	a_hgvs_c_column: str,
	b_hgvs_c_column: str,
	a_hgvs_g_column: str,
	b_hgvs_g_column: str,
	compare: str,
	match_columns: tuple[str, ...],
	ignore_reference_differences: bool,
	a_prefix: str,
	b_prefix: str,
	pass_through_selected_columns: tuple[str, ...],
	pass_through_all_columns: bool,
	a_pass_through_prefix: str,
	b_pass_through_prefix: str,
	a_only_output: Path,
	b_only_output: Path,
	different_output: Path,
	a_missing_output: Path,
	b_missing_output: Path,
) -> None:
	a_input_format = a_input_format.lower()
	b_input_format = b_input_format.lower()
	compare = compare.lower()
	configure_csv_field_size_limit()
	match_column_names = parse_column_list_option(match_columns, "--match-columns")
	selected_pass_through_column_names = parse_column_list_option(
		pass_through_selected_columns,
		"--pass-through-columns",
	)
	if pass_through_all_columns and selected_pass_through_column_names:
		raise click.ClickException(
			"--pass-through-all-columns cannot be combined with --pass-through-columns. Use one strategy."
		)

	for option_name, option_value in (
		("--a-separator", a_separator),
		("--b-separator", b_separator),
		("--a-prefix", a_prefix),
		("--b-prefix", b_prefix),
		("--a-pass-through-prefix", a_pass_through_prefix),
		("--b-pass-through-prefix", b_pass_through_prefix),
	):
		if not option_value:
			raise click.ClickException(f"{option_name} cannot be empty.")
		if "\t" in option_value:
			raise click.ClickException(f"{option_name} cannot contain a tab character.")

	a_spec = InputSpec(
		label="A",
		path=a_input,
		input_format=a_input_format,
		separator=a_separator,
		transcript_column=a_transcript_column,
		hgvs_p_column=a_hgvs_p_column,
		hgvs_c_column=a_hgvs_c_column,
		hgvs_g_column=a_hgvs_g_column,
	)
	b_spec = InputSpec(
		label="B",
		path=b_input,
		input_format=b_input_format,
		separator=b_separator,
		transcript_column=b_transcript_column,
		hgvs_p_column=b_hgvs_p_column,
		hgvs_c_column=b_hgvs_c_column,
		hgvs_g_column=b_hgvs_g_column,
	)

	a_field_names, a_rows_by_key, a_key_order, a_missing_key_rows = read_rows_by_key(
		a_spec,
		compare=compare,
		match_columns=match_column_names,
	)
	b_field_names, b_rows_by_key, b_key_order, b_missing_key_rows = read_rows_by_key(
		b_spec,
		compare=compare,
		match_columns=match_column_names,
	)

	a_keys = set(a_rows_by_key)
	b_keys = set(b_rows_by_key)

	a_only_rows = [a_rows_by_key[key] for key in a_key_order if key not in b_keys]
	b_only_rows = [b_rows_by_key[key] for key in b_key_order if key not in a_keys]

	a_core_columns = [
		column
		for column in [a_spec.transcript_column, a_spec.hgvs_p_column, a_spec.hgvs_c_column, a_spec.hgvs_g_column]
		if column in a_field_names
	]
	b_core_columns = [
		column
		for column in [b_spec.transcript_column, b_spec.hgvs_p_column, b_spec.hgvs_c_column, b_spec.hgvs_g_column]
		if column in b_field_names
	]

	a_column_map = build_prefixed_column_names(a_core_columns, a_prefix)
	b_column_map = build_prefixed_column_names(b_core_columns, b_prefix)

	a_additional_columns = [column for column in a_field_names if column not in a_core_columns]
	b_additional_columns = [column for column in b_field_names if column not in b_core_columns]

	a_selected_additional_columns: list[str] = []
	b_selected_additional_columns: list[str] = []

	if pass_through_all_columns:
		a_selected_additional_columns = list(a_additional_columns)
		b_selected_additional_columns = list(b_additional_columns)
	else:
		unknown_pass_through_columns: list[str] = []
		for column_name in selected_pass_through_column_names:
			found_in_at_least_one_input = False
			if column_name in a_additional_columns:
				a_selected_additional_columns.append(column_name)
				found_in_at_least_one_input = True
			if column_name in b_additional_columns:
				b_selected_additional_columns.append(column_name)
				found_in_at_least_one_input = True
			if not found_in_at_least_one_input:
				unknown_pass_through_columns.append(column_name)

		if unknown_pass_through_columns:
			raise click.ClickException(
				"Requested --pass-through-columns not found among non-core columns in A or B: "
				+ ", ".join(unknown_pass_through_columns)
			)

	a_extra_column_map: dict[str, str] = {}
	b_extra_column_map: dict[str, str] = {}
	if pass_through_all_columns or selected_pass_through_column_names:
		a_extra_column_map = build_prefixed_column_names(a_selected_additional_columns, a_pass_through_prefix)
		b_extra_column_map = build_prefixed_column_names(b_selected_additional_columns, b_pass_through_prefix)

	differences_field_names = list(a_column_map.values()) + list(b_column_map.values())
	if pass_through_all_columns or selected_pass_through_column_names:
		differences_field_names.extend(a_extra_column_map.values())
		differences_field_names.extend(b_extra_column_map.values())

	if len(set(differences_field_names)) != len(differences_field_names):
		raise click.ClickException(
			"Output column name collision in differences output. "
			"Adjust --a-prefix/--b-prefix and/or pass-through prefixes."
		)

	difference_rows: list[dict[str, str]] = []
	for row_key in a_key_order:
		if row_key not in b_rows_by_key:
			continue

		a_row = a_rows_by_key[row_key]
		b_row = b_rows_by_key[row_key]
		if not rows_differ(
			a_row,
			b_row,
			a_spec=a_spec,
			b_spec=b_spec,
			compare=compare,
			ignore_reference_differences=ignore_reference_differences,
		):
			continue

		difference_row = {
			**{output_name: a_row.get(input_name, "") for input_name, output_name in a_column_map.items()},
			**{output_name: b_row.get(input_name, "") for input_name, output_name in b_column_map.items()},
		}
		if pass_through_all_columns or selected_pass_through_column_names:
			difference_row.update(
				{output_name: a_row.get(input_name, "") for input_name, output_name in a_extra_column_map.items()}
			)
			difference_row.update(
				{output_name: b_row.get(input_name, "") for input_name, output_name in b_extra_column_map.items()}
			)
		difference_rows.append(difference_row)

	write_tsv(a_only_output, field_names=a_field_names, rows=a_only_rows)
	write_tsv(b_only_output, field_names=b_field_names, rows=b_only_rows)
	write_tsv(different_output, field_names=differences_field_names, rows=difference_rows)
	write_tsv(a_missing_output, field_names=a_field_names, rows=a_missing_key_rows)
	write_tsv(b_missing_output, field_names=b_field_names, rows=b_missing_key_rows)

	click.echo(
		"Comparison summary: "
		f"A rows={len(a_rows_by_key)}, "
		f"B rows={len(b_rows_by_key)}, "
		f"A missing transcript/hgvs_p={len(a_missing_key_rows)}, "
		f"B missing transcript/hgvs_p={len(b_missing_key_rows)}, "
		f"A-only={len(a_only_rows)}, "
		f"B-only={len(b_only_rows)}, "
		f"different={len(difference_rows)}.",
		err=True,
	)


if __name__ == "__main__":
	main()
