import csv
from pathlib import Path
from unittest.mock import Mock, patch

from click.testing import CliRunner
import pytest

from src.scripts import reverse_translate_variants as rtv


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=",")
        writer.writeheader()
        writer.writerows(rows)


def test_parse_hgvs_protein_change() -> None:
    parsed = rtv.parse_hgvs_protein_change("p.Arg175His")
    assert parsed.reference_aa == "R"
    assert parsed.position == 175
    assert parsed.alternate_aa == "H"

    parsed_synonymous = rtv.parse_hgvs_protein_change("NP_000537.3:p.Arg175=")
    assert parsed_synonymous.reference_aa == "R"
    assert parsed_synonymous.position == 175
    assert parsed_synonymous.alternate_aa == "R"


def test_reverse_translate_hgvs_p_with_mocked_data_provider() -> None:
    data_provider = Mock()
    data_provider.get_tx_identity_info.return_value = {"cds_start_i": 0, "cds_end_i": 3}
    data_provider.get_seq.return_value = "ATG"

    transcript_cache: dict[str, tuple[str, str]] = {}

    with patch(
        "src.scripts.reverse_translate_variants.map_hgvs_c_to_hgvs_g",
        side_effect=lambda parser, mapper, tx, hgvs_c: f"GENOMIC:{tx}:{hgvs_c}",
    ):
        rows = rtv.reverse_translate_hgvs_p(
            transcript_accession="NM_TEST.1",
            hgvs_protein="p.Met1Val",
            include_indels=False,
            max_indel_size=3,
            strict_ref_aa=True,
            use_inv_notation=False,
            allow_length_changing_stop_candidates=True,
            parser=Mock(),
            mapper=Mock(),
            data_provider=data_provider,
            transcript_cache=transcript_cache,
        )

    assert len(rows) == 1
    assert rows[0]["variant_type"] == "snv"
    assert rows[0]["hgvs_c"] == "NM_TEST.1:c.1A>G"
    assert rows[0]["hgvs_g"] == "GENOMIC:NM_TEST.1:c.1A>G"


def test_reverse_translate_non_snv_substitution_includes_triplet_delins() -> None:
    data_provider = Mock()
    data_provider.get_tx_identity_info.return_value = {"cds_start_i": 0, "cds_end_i": 3}
    data_provider.get_seq.return_value = "GCT"  # Ala codon

    transcript_cache: dict[str, tuple[str, str]] = {}

    with patch(
        "src.scripts.reverse_translate_variants.map_hgvs_c_to_hgvs_g",
        side_effect=lambda parser, mapper, tx, hgvs_c: f"GENOMIC:{tx}:{hgvs_c}",
    ):
        rows = rtv.reverse_translate_hgvs_p(
            transcript_accession="NM_TEST.1",
            hgvs_protein="p.Ala1Trp",  # Requires multi-nt change; not SNV-accessible from GCT
            include_indels=True,
            max_indel_size=3,
            strict_ref_aa=True,
            use_inv_notation=False,
            allow_length_changing_stop_candidates=True,
            parser=Mock(),
            mapper=Mock(),
            data_provider=data_provider,
            transcript_cache=transcript_cache,
        )

    assert rows
    assert any(row["variant_type"] == "delins" for row in rows)
    assert any(row["hgvs_c"] == "NM_TEST.1:c.1_3delinsTGG" for row in rows)
    assert not any(row["hgvs_c"] == "NM_TEST.1:c.1_2delinsTGG" for row in rows)
    assert not any(row["hgvs_c"] == "NM_TEST.1:c.1delinsTGG" for row in rows)


def test_reverse_translate_prefers_minimal_two_nt_delins_over_full_codon_replacement() -> None:
    data_provider = Mock()
    data_provider.get_tx_identity_info.return_value = {"cds_start_i": 0, "cds_end_i": 3}
    data_provider.get_seq.return_value = "GCT"  # Ala codon

    transcript_cache: dict[str, tuple[str, str]] = {}

    with patch(
        "src.scripts.reverse_translate_variants.map_hgvs_c_to_hgvs_g",
        side_effect=lambda parser, mapper, tx, hgvs_c: f"GENOMIC:{tx}:{hgvs_c}",
    ):
        rows = rtv.reverse_translate_hgvs_p(
            transcript_accession="NM_TEST.1",
            hgvs_protein="p.Ala1Leu",  # Reachable via 2-nt delins such as GCT -> CTT
            include_indels=True,
            max_indel_size=3,
            strict_ref_aa=True,
            use_inv_notation=False,
            allow_length_changing_stop_candidates=True,
            parser=Mock(),
            mapper=Mock(),
            data_provider=data_provider,
            transcript_cache=transcript_cache,
        )

    assert rows
    assert any(row["hgvs_c"] == "NM_TEST.1:c.1_2delinsCT" for row in rows)
    assert not any(row["hgvs_c"] == "NM_TEST.1:c.1_3delinsCTT" for row in rows)
    assert not any(row["hgvs_c"] == "NM_TEST.1:c.1delinsCT" for row in rows)


def test_reverse_translate_can_emit_inv_for_full_codon_inversion() -> None:
    data_provider = Mock()
    data_provider.get_tx_identity_info.return_value = {"cds_start_i": 0, "cds_end_i": 3}
    data_provider.get_seq.return_value = "GCT"  # Ala codon; reverse complement is AGC (Ser)

    transcript_cache: dict[str, tuple[str, str]] = {}

    with patch(
        "src.scripts.reverse_translate_variants.map_hgvs_c_to_hgvs_g",
        side_effect=lambda parser, mapper, tx, hgvs_c: f"GENOMIC:{tx}:{hgvs_c}",
    ):
        rows = rtv.reverse_translate_hgvs_p(
            transcript_accession="NM_TEST.1",
            hgvs_protein="p.Ala1Ser",
            include_indels=True,
            max_indel_size=3,
            strict_ref_aa=True,
            use_inv_notation=True,
            allow_length_changing_stop_candidates=True,
            parser=Mock(),
            mapper=Mock(),
            data_provider=data_provider,
            transcript_cache=transcript_cache,
        )

    assert any(row["variant_type"] == "inv" and row["hgvs_c"] == "NM_TEST.1:c.1_3inv" for row in rows)
    assert not any(row["hgvs_c"] == "NM_TEST.1:c.1_3delinsAGC" for row in rows)


def test_reverse_translate_can_emit_inv_for_minimized_two_nt_inversion() -> None:
    data_provider = Mock()
    data_provider.get_tx_identity_info.return_value = {"cds_start_i": 0, "cds_end_i": 3}
    data_provider.get_seq.return_value = "GCT"  # Positions 2-3: CT -> AG (reverse complement), yielding GAG (Glu)

    transcript_cache: dict[str, tuple[str, str]] = {}

    with patch(
        "src.scripts.reverse_translate_variants.map_hgvs_c_to_hgvs_g",
        side_effect=lambda parser, mapper, tx, hgvs_c: f"GENOMIC:{tx}:{hgvs_c}",
    ):
        rows = rtv.reverse_translate_hgvs_p(
            transcript_accession="NM_TEST.1",
            hgvs_protein="p.Ala1Glu",
            include_indels=True,
            max_indel_size=3,
            strict_ref_aa=True,
            use_inv_notation=True,
            allow_length_changing_stop_candidates=True,
            parser=Mock(),
            mapper=Mock(),
            data_provider=data_provider,
            transcript_cache=transcript_cache,
        )

    assert any(row["variant_type"] == "inv" and row["hgvs_c"] == "NM_TEST.1:c.2_3inv" for row in rows)
    assert not any(row["hgvs_c"] == "NM_TEST.1:c.2_3delinsAG" for row in rows)


def test_reverse_translate_premature_stop_includes_length_changing_candidates_by_default() -> None:
    insertion_data_provider = Mock()
    insertion_data_provider.get_tx_identity_info.return_value = {"cds_start_i": 0, "cds_end_i": 6}
    insertion_data_provider.get_seq.return_value = "TTTGCT"  # Phe Ala

    deletion_data_provider = Mock()
    deletion_data_provider.get_tx_identity_info.return_value = {"cds_start_i": 0, "cds_end_i": 6}
    deletion_data_provider.get_seq.return_value = "TTTTAA"  # Phe Stop

    with patch(
        "src.scripts.reverse_translate_variants.map_hgvs_c_to_hgvs_g",
        side_effect=lambda parser, mapper, tx, hgvs_c: f"GENOMIC:{tx}:{hgvs_c}",
    ):
        insertion_rows = rtv.reverse_translate_hgvs_p(
            transcript_accession="NM_INS.1",
            hgvs_protein="p.Phe1Ter",
            include_indels=True,
            max_indel_size=3,
            strict_ref_aa=True,
            use_inv_notation=False,
            allow_length_changing_stop_candidates=True,
            parser=Mock(),
            mapper=Mock(),
            data_provider=insertion_data_provider,
            transcript_cache={},
        )
        deletion_rows = rtv.reverse_translate_hgvs_p(
            transcript_accession="NM_DEL.1",
            hgvs_protein="p.Phe1Ter",
            include_indels=True,
            max_indel_size=3,
            strict_ref_aa=True,
            use_inv_notation=False,
            allow_length_changing_stop_candidates=True,
            parser=Mock(),
            mapper=Mock(),
            data_provider=deletion_data_provider,
            transcript_cache={},
        )

    assert any(row["variant_type"] == "insertion" and row["hgvs_c"] == "NM_INS.1:c.1_2insAAA" for row in insertion_rows)
    assert any(row["variant_type"] == "deletion" and row["hgvs_c"] == "NM_DEL.1:c.1_3del" for row in deletion_rows)


def test_reverse_translate_premature_stop_can_restrict_to_substitutions_and_same_length_delins() -> None:
    insertion_data_provider = Mock()
    insertion_data_provider.get_tx_identity_info.return_value = {"cds_start_i": 0, "cds_end_i": 6}
    insertion_data_provider.get_seq.return_value = "TTTGCT"  # Phe Ala

    deletion_data_provider = Mock()
    deletion_data_provider.get_tx_identity_info.return_value = {"cds_start_i": 0, "cds_end_i": 6}
    deletion_data_provider.get_seq.return_value = "TTTTAA"  # Phe Stop

    with patch(
        "src.scripts.reverse_translate_variants.map_hgvs_c_to_hgvs_g",
        side_effect=lambda parser, mapper, tx, hgvs_c: f"GENOMIC:{tx}:{hgvs_c}",
    ):
        insertion_rows = rtv.reverse_translate_hgvs_p(
            transcript_accession="NM_INS.1",
            hgvs_protein="p.Phe1Ter",
            include_indels=True,
            max_indel_size=3,
            strict_ref_aa=True,
            use_inv_notation=False,
            allow_length_changing_stop_candidates=False,
            parser=Mock(),
            mapper=Mock(),
            data_provider=insertion_data_provider,
            transcript_cache={},
        )
        deletion_rows = rtv.reverse_translate_hgvs_p(
            transcript_accession="NM_DEL.1",
            hgvs_protein="p.Phe1Ter",
            include_indels=True,
            max_indel_size=3,
            strict_ref_aa=True,
            use_inv_notation=False,
            allow_length_changing_stop_candidates=False,
            parser=Mock(),
            mapper=Mock(),
            data_provider=deletion_data_provider,
            transcript_cache={},
        )

    assert insertion_rows
    assert any(row["variant_type"] == "delins" and row["hgvs_c"] == "NM_INS.1:c.2_3delinsAA" for row in insertion_rows)
    assert not any(row["variant_type"] == "insertion" for row in insertion_rows)
    assert not any(row["variant_type"] == "deletion" for row in insertion_rows)
    assert not any(row["variant_type"] == "insertion" for row in deletion_rows)
    assert not any(row["variant_type"] == "deletion" for row in deletion_rows)


def test_join_variant_rows_keeps_aligned_counts() -> None:
    joined = rtv.join_variant_rows(
        [
            {"variant_type": "snv", "hgvs_c": "NM_TEST.1:c.1A>G", "hgvs_g": "NC_000001.11:g.1A>G"},
            {"variant_type": "delins", "hgvs_c": "NM_TEST.1:c.1_3delinsTGG", "hgvs_g": "NC_000001.11:g.1_3delinsTGG"},
        ],
        join_delimiter="|",
    )

    assert joined["variant_type"] == "snv|delins"
    assert joined["hgvs_c"] == "NM_TEST.1:c.1A>G|NM_TEST.1:c.1_3delinsTGG"
    assert joined["hgvs_g"] == "NC_000001.11:g.1A>G|NC_000001.11:g.1_3delinsTGG"
    assert len(joined["hgvs_c"].split("|")) == len(joined["hgvs_g"].split("|"))


def test_reverse_translate_batch_rows_one_row_per_input_programmatic_helper() -> None:
    rows = [
        {"id": "r1", "transcript": "NM_A.1", "hgvs_p": "R2H"},
        {"id": "r2", "transcript": "NM_B.1", "hgvs_p": ""},
        {"id": "r3", "transcript": "", "hgvs_p": "p.Gly2="},
    ]

    reverse_translate_calls: list[str] = []

    def mocked_reverse_translate(**kwargs):
        reverse_translate_calls.append(kwargs["hgvs_protein"])
        return [
            {
                "variant_type": "snv",
                "hgvs_c": "NM_A.1:c.10A>G",
                "hgvs_g": "NC_000001.11:g.100A>G",
            },
            {
                "variant_type": "delins",
                "hgvs_c": "NM_A.1:c.10_11delinsTT",
                "hgvs_g": "NC_000001.11:g.100_101delinsTT",
            },
        ]

    with patch("src.scripts.reverse_translate_variants.reverse_translate_hgvs_p", side_effect=mocked_reverse_translate):
        output_rows = rtv.reverse_translate_batch_rows(
            rows=rows,
            transcript_column="transcript",
            hgvs_p_column="hgvs_p",
            include_indels=True,
            max_indel_size=3,
            strict_ref_aa=True,
            use_inv_notation=False,
            allow_length_changing_stop_candidates=True,
            parser=Mock(),
            mapper=Mock(),
            data_provider=Mock(),
            transcript_cache={},
            auto_format_hgvs_p=True,
            one_row_per_input=True,
            join_delimiter="^",
        )

    assert reverse_translate_calls == ["p.R2H"]
    assert len(output_rows) == 2
    assert output_rows[0]["id"] == "r1"
    assert output_rows[0]["variant_type"] == "snv^delins"
    assert output_rows[0]["hgvs_c"] == "NM_A.1:c.10A>G^NM_A.1:c.10_11delinsTT"
    assert output_rows[1]["id"] == "r3"
    assert output_rows[1]["variant_type"] == ""


def test_reverse_translate_batch_rows_can_keep_missing_hgvs_p_rows() -> None:
    rows = [{"id": "r1", "transcript": "NM_A.1", "hgvs_p": ""}]

    output_rows = rtv.reverse_translate_batch_rows(
        rows=rows,
        transcript_column="transcript",
        hgvs_p_column="hgvs_p",
        include_indels=False,
        max_indel_size=3,
        strict_ref_aa=True,
        use_inv_notation=False,
        allow_length_changing_stop_candidates=True,
        parser=Mock(),
        mapper=Mock(),
        data_provider=Mock(),
        transcript_cache={},
        skip_missing_hgvs_p=False,
    )

    assert len(output_rows) == 1
    assert output_rows[0]["id"] == "r1"
    assert output_rows[0]["variant_type"] == ""
    assert output_rows[0]["hgvs_c"] == ""
    assert output_rows[0]["hgvs_g"] == ""


def test_reverse_translate_batch_rows_can_raise_on_error() -> None:
    rows = [{"id": "r1", "transcript": "NM_A.1", "hgvs_p": "p.Arg2His"}]

    with patch(
        "src.scripts.reverse_translate_variants.reverse_translate_hgvs_p",
        side_effect=RuntimeError("boom"),
    ):
        with pytest.raises(RuntimeError, match="boom"):
            rtv.reverse_translate_batch_rows(
                rows=rows,
                transcript_column="transcript",
                hgvs_p_column="hgvs_p",
                include_indels=False,
                max_indel_size=3,
                strict_ref_aa=True,
                use_inv_notation=False,
                allow_length_changing_stop_candidates=True,
                parser=Mock(),
                mapper=Mock(),
                data_provider=Mock(),
                transcript_cache={},
                raise_on_error=True,
            )


def test_reverse_translate_batch_rows_can_collect_error_rows() -> None:
    rows = [
        {"id": "r1", "transcript": "", "hgvs_p": "p.Arg2His"},
        {"id": "r2", "transcript": "NM_A.1", "hgvs_p": "p.Arg2His"},
    ]
    error_rows: list[dict[str, str]] = []

    with patch(
        "src.scripts.reverse_translate_variants.reverse_translate_hgvs_p",
        side_effect=RuntimeError("boom"),
    ):
        output_rows = rtv.reverse_translate_batch_rows(
            rows=rows,
            transcript_column="transcript",
            hgvs_p_column="hgvs_p",
            include_indels=False,
            max_indel_size=3,
            strict_ref_aa=True,
            use_inv_notation=False,
            allow_length_changing_stop_candidates=True,
            parser=Mock(),
            mapper=Mock(),
            data_provider=Mock(),
            transcript_cache={},
            error_rows=error_rows,
        )

    assert len(output_rows) == 2
    assert len(error_rows) == 2
    assert error_rows[0]["id"] == "r1"
    assert error_rows[0]["error"] == "Missing transcript value."
    assert error_rows[1]["id"] == "r2"
    assert error_rows[1]["error"] == "Failed reverse translation (boom)"


def test_batch_mode_pass_through_all_columns_with_prefix(runner: CliRunner, tmp_path: Path) -> None:
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"

    write_tsv(
        input_path,
        rows=[
            {
                "sample_id": "row1",
                "note": "keep me",
                "transcript": "NM_A.1",
                "hgvs_p": "p.Arg2His",
            },
            {
                "sample_id": "row2",
                "note": "no variants",
                "transcript": "NM_B.1",
                "hgvs_p": "p.Gly2=",
            },
        ],
        fieldnames=["sample_id", "note", "transcript", "hgvs_p"],
    )

    mocked_connect = Mock()
    mocked_mapper = Mock()

    def mocked_reverse_translate(**kwargs):
        if kwargs["transcript_accession"] == "NM_A.1":
            return [
                {
                    "variant_type": "snv",
                    "hgvs_c": "NM_A.1:c.10A>G",
                    "hgvs_g": "NC_000001.11:g.100A>G",
                },
                {
                    "variant_type": "snv",
                    "hgvs_c": "NM_A.1:c.10A>T",
                    "hgvs_g": "NC_000001.11:g.100A>T",
                },
            ]
        return []

    with (
        patch("src.scripts.reverse_translate_variants.hgvs.dataproviders.uta.connect", return_value=mocked_connect),
        patch("src.scripts.reverse_translate_variants.hgvs.assemblymapper.AssemblyMapper", return_value=mocked_mapper),
        patch("src.scripts.reverse_translate_variants.reverse_translate_hgvs_p", side_effect=mocked_reverse_translate),
    ):
        result = runner.invoke(
            rtv.main,
            [
                "--input",
                str(input_path),
                "--pass-through-all-columns",
                "--pass-through-prefix",
                "meta_",
                "--output",
                str(output_path),
            ],
        )

    assert result.exit_code == 0, result.output

    with output_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    assert len(rows) == 3
    assert reader.fieldnames is not None
    assert "meta_sample_id" in reader.fieldnames
    assert "meta_note" in reader.fieldnames
    assert "variant_type" in reader.fieldnames
    assert "hgvs_c" in reader.fieldnames
    assert "hgvs_g" in reader.fieldnames

    assert rows[0]["meta_sample_id"] == "row1"
    assert rows[0]["meta_note"] == "keep me"
    assert rows[0]["hgvs_c"] == "NM_A.1:c.10A>G"

    assert rows[1]["meta_sample_id"] == "row1"
    assert rows[1]["meta_note"] == "keep me"
    assert rows[1]["hgvs_c"] == "NM_A.1:c.10A>T"

    assert rows[2]["meta_sample_id"] == "row2"
    assert rows[2]["meta_note"] == "no variants"
    assert rows[2]["variant_type"] == ""
    assert rows[2]["hgvs_c"] == ""
    assert rows[2]["hgvs_g"] == ""


def test_batch_mode_pass_through_selected_columns_only(runner: CliRunner, tmp_path: Path) -> None:
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"

    write_tsv(
        input_path,
        rows=[
            {
                "sample_id": "row1",
                "note": "keep me",
                "batch": "b1",
                "transcript": "NM_A.1",
                "hgvs_p": "p.Arg2His",
            }
        ],
        fieldnames=["sample_id", "note", "batch", "transcript", "hgvs_p"],
    )

    with (
        patch("src.scripts.reverse_translate_variants.hgvs.dataproviders.uta.connect", return_value=Mock()),
        patch("src.scripts.reverse_translate_variants.hgvs.assemblymapper.AssemblyMapper", return_value=Mock()),
        patch(
            "src.scripts.reverse_translate_variants.reverse_translate_hgvs_p",
            return_value=[
                {
                    "variant_type": "snv",
                    "hgvs_c": "NM_A.1:c.10A>G",
                    "hgvs_g": "NC_000001.11:g.100A>G",
                }
            ],
        ),
    ):
        result = runner.invoke(
            rtv.main,
            [
                "--input",
                str(input_path),
                "--pass-through-columns",
                "sample_id",
                "--pass-through-columns",
                "batch",
                "--pass-through-prefix",
                "meta_",
                "--output",
                str(output_path),
            ],
        )

    assert result.exit_code == 0, result.output

    with output_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    assert len(rows) == 1
    assert reader.fieldnames is not None
    assert "meta_sample_id" in reader.fieldnames
    assert "meta_batch" in reader.fieldnames
    assert "meta_note" not in reader.fieldnames
    assert rows[0]["meta_sample_id"] == "row1"
    assert rows[0]["meta_batch"] == "b1"


def test_batch_mode_pass_through_columns_accepts_csv_list_and_repeated_uses(
    runner: CliRunner,
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"

    write_tsv(
        input_path,
        rows=[
            {
                "sample_id": "row1",
                "note": "keep me",
                "batch": "b1",
                "plate": "p7",
                "transcript": "NM_A.1",
                "hgvs_p": "p.Arg2His",
            }
        ],
        fieldnames=["sample_id", "note", "batch", "plate", "transcript", "hgvs_p"],
    )

    with (
        patch("src.scripts.reverse_translate_variants.hgvs.dataproviders.uta.connect", return_value=Mock()),
        patch("src.scripts.reverse_translate_variants.hgvs.assemblymapper.AssemblyMapper", return_value=Mock()),
        patch(
            "src.scripts.reverse_translate_variants.reverse_translate_hgvs_p",
            return_value=[
                {
                    "variant_type": "snv",
                    "hgvs_c": "NM_A.1:c.10A>G",
                    "hgvs_g": "NC_000001.11:g.100A>G",
                }
            ],
        ),
    ):
        result = runner.invoke(
            rtv.main,
            [
                "--input",
                str(input_path),
                "--pass-through-columns",
                "sample_id,batch",
                "--pass-through-columns",
                "plate",
                "--pass-through-prefix",
                "meta_",
                "--output",
                str(output_path),
            ],
        )

    assert result.exit_code == 0, result.output

    with output_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    assert len(rows) == 1
    assert reader.fieldnames is not None
    assert "meta_sample_id" in reader.fieldnames
    assert "meta_batch" in reader.fieldnames
    assert "meta_plate" in reader.fieldnames
    assert "meta_note" not in reader.fieldnames


def test_batch_mode_does_not_pass_through_additional_columns_by_default(runner: CliRunner, tmp_path: Path) -> None:
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"

    write_tsv(
        input_path,
        rows=[
            {
                "sample_id": "row1",
                "note": "keep me",
                "transcript": "NM_A.1",
                "hgvs_p": "p.Arg2His",
            }
        ],
        fieldnames=["sample_id", "note", "transcript", "hgvs_p"],
    )

    with (
        patch("src.scripts.reverse_translate_variants.hgvs.dataproviders.uta.connect", return_value=Mock()),
        patch("src.scripts.reverse_translate_variants.hgvs.assemblymapper.AssemblyMapper", return_value=Mock()),
        patch(
            "src.scripts.reverse_translate_variants.reverse_translate_hgvs_p",
            return_value=[
                {
                    "variant_type": "snv",
                    "hgvs_c": "NM_A.1:c.10A>G",
                    "hgvs_g": "NC_000001.11:g.100A>G",
                }
            ],
        ),
    ):
        result = runner.invoke(
            rtv.main,
            [
                "--input",
                str(input_path),
                "--output",
                str(output_path),
            ],
        )

    assert result.exit_code == 0, result.output

    with output_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    assert len(rows) == 1
    assert reader.fieldnames is not None
    assert "sample_id" not in reader.fieldnames
    assert "note" not in reader.fieldnames
    assert "transcript" in reader.fieldnames
    assert "hgvs_p" in reader.fieldnames
    assert rows[0]["transcript"] == "NM_A.1"
    assert rows[0]["hgvs_c"] == "NM_A.1:c.10A>G"


def test_batch_mode_fails_when_transcript_column_missing(runner: CliRunner, tmp_path: Path) -> None:
    input_path = tmp_path / "missing_transcript.tsv"

    write_tsv(
        input_path,
        rows=[{"sample_id": "row1", "hgvs_p": "p.Arg175His"}],
        fieldnames=["sample_id", "hgvs_p"],
    )

    with (
        patch("src.scripts.reverse_translate_variants.hgvs.dataproviders.uta.connect", return_value=Mock()),
        patch("src.scripts.reverse_translate_variants.hgvs.assemblymapper.AssemblyMapper", return_value=Mock()),
    ):
        result = runner.invoke(
            rtv.main,
            [
                "--input",
                str(input_path),
            ],
        )

    assert result.exit_code != 0
    assert "Missing transcript column in input TSV" in result.output


def test_batch_mode_fails_when_hgvs_p_column_missing(runner: CliRunner, tmp_path: Path) -> None:
    input_path = tmp_path / "missing_hgvs_p.tsv"

    write_tsv(
        input_path,
        rows=[{"sample_id": "row1", "transcript": "NM_000546.6"}],
        fieldnames=["sample_id", "transcript"],
    )

    with (
        patch("src.scripts.reverse_translate_variants.hgvs.dataproviders.uta.connect", return_value=Mock()),
        patch("src.scripts.reverse_translate_variants.hgvs.assemblymapper.AssemblyMapper", return_value=Mock()),
    ):
        result = runner.invoke(
            rtv.main,
            [
                "--input",
                str(input_path),
            ],
        )

    assert result.exit_code != 0
    assert "Missing HGVS p column in input TSV" in result.output


def test_batch_mode_skips_rows_missing_hgvs_p_with_warning(runner: CliRunner, tmp_path: Path) -> None:
    input_path = tmp_path / "missing_hgvs_p_values.tsv"
    output_path = tmp_path / "output.tsv"

    write_tsv(
        input_path,
        rows=[
            {"transcript": "NM_A.1", "hgvs_p": "p.Arg2His"},
            {"transcript": "NM_B.1", "hgvs_p": ""},
            {"transcript": "NM_C.1", "hgvs_p": "p.Gly2="},
        ],
        fieldnames=["transcript", "hgvs_p"],
    )

    mocked_connect = Mock()
    mocked_mapper = Mock()
    mocked_reverse_translate = Mock(
        return_value=[
            {
                "variant_type": "snv",
                "hgvs_c": "NM_TEST.1:c.10A>G",
                "hgvs_g": "NC_000001.11:g.100A>G",
            }
        ]
    )

    with (
        patch("src.scripts.reverse_translate_variants.hgvs.dataproviders.uta.connect", return_value=mocked_connect),
        patch("src.scripts.reverse_translate_variants.hgvs.assemblymapper.AssemblyMapper", return_value=mocked_mapper),
        patch("src.scripts.reverse_translate_variants.reverse_translate_hgvs_p", mocked_reverse_translate),
    ):
        result = runner.invoke(
            rtv.main,
            [
                "--input",
                str(input_path),
                "--output",
                str(output_path),
            ],
        )

    assert result.exit_code == 0, result.output
    assert mocked_reverse_translate.call_count == 2
    assert "Warning: Skipped 1 row(s) with missing HGVS p. values." in result.output

    with output_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    assert len(rows) == 2
    assert [row["transcript"] for row in rows] == ["NM_A.1", "NM_C.1"]


def test_batch_mode_limit_processes_first_n_rows(runner: CliRunner, tmp_path: Path) -> None:
    input_path = tmp_path / "limit_input.tsv"
    output_path = tmp_path / "limit_output.tsv"

    write_tsv(
        input_path,
        rows=[
            {"sample_id": "row1", "transcript": "NM_A.1", "hgvs_p": "p.Arg2His"},
            {"sample_id": "row2", "transcript": "NM_B.1", "hgvs_p": "p.Arg2His"},
            {"sample_id": "row3", "transcript": "NM_C.1", "hgvs_p": "p.Arg2His"},
        ],
        fieldnames=["sample_id", "transcript", "hgvs_p"],
    )

    mocked_connect = Mock()
    mocked_mapper = Mock()

    def mocked_reverse_translate(**kwargs):
        return [
            {
                "variant_type": "snv",
                "hgvs_c": f"{kwargs['transcript_accession']}:c.10A>G",
                "hgvs_g": "NC_000001.11:g.100A>G",
            }
        ]

    with (
        patch("src.scripts.reverse_translate_variants.hgvs.dataproviders.uta.connect", return_value=mocked_connect),
        patch("src.scripts.reverse_translate_variants.hgvs.assemblymapper.AssemblyMapper", return_value=mocked_mapper),
        patch("src.scripts.reverse_translate_variants.reverse_translate_hgvs_p", side_effect=mocked_reverse_translate),
    ):
        result = runner.invoke(
            rtv.main,
            [
                "--input",
                str(input_path),
                "--limit",
                "2",
                "--output",
                str(output_path),
            ],
        )

    assert result.exit_code == 0, result.output

    with output_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    assert len(rows) == 2
    assert [row["transcript"] for row in rows] == ["NM_A.1", "NM_B.1"]


def test_batch_mode_skip_skips_first_n_rows(runner: CliRunner, tmp_path: Path) -> None:
    input_path = tmp_path / "skip_input.tsv"
    output_path = tmp_path / "skip_output.tsv"

    write_tsv(
        input_path,
        rows=[
            {"sample_id": "row1", "transcript": "NM_A.1", "hgvs_p": "p.Arg2His"},
            {"sample_id": "row2", "transcript": "NM_B.1", "hgvs_p": "p.Arg2His"},
            {"sample_id": "row3", "transcript": "NM_C.1", "hgvs_p": "p.Arg2His"},
        ],
        fieldnames=["sample_id", "transcript", "hgvs_p"],
    )

    mocked_connect = Mock()
    mocked_mapper = Mock()

    def mocked_reverse_translate(**kwargs):
        return [
            {
                "variant_type": "snv",
                "hgvs_c": f"{kwargs['transcript_accession']}:c.10A>G",
                "hgvs_g": "NC_000001.11:g.100A>G",
            }
        ]

    with (
        patch("src.scripts.reverse_translate_variants.hgvs.dataproviders.uta.connect", return_value=mocked_connect),
        patch("src.scripts.reverse_translate_variants.hgvs.assemblymapper.AssemblyMapper", return_value=mocked_mapper),
        patch("src.scripts.reverse_translate_variants.reverse_translate_hgvs_p", side_effect=mocked_reverse_translate),
    ):
        result = runner.invoke(
            rtv.main,
            [
                "--input",
                str(input_path),
                "--skip",
                "1",
                "--output",
                str(output_path),
            ],
        )

    assert result.exit_code == 0, result.output

    with output_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    assert len(rows) == 2
    assert [row["transcript"] for row in rows] == ["NM_B.1", "NM_C.1"]


def test_batch_mode_skip_is_applied_before_limit(runner: CliRunner, tmp_path: Path) -> None:
    input_path = tmp_path / "skip_limit_input.tsv"
    output_path = tmp_path / "skip_limit_output.tsv"

    write_tsv(
        input_path,
        rows=[
            {"sample_id": "row1", "transcript": "NM_A.1", "hgvs_p": "p.Arg2His"},
            {"sample_id": "row2", "transcript": "NM_B.1", "hgvs_p": "p.Arg2His"},
            {"sample_id": "row3", "transcript": "NM_C.1", "hgvs_p": "p.Arg2His"},
        ],
        fieldnames=["sample_id", "transcript", "hgvs_p"],
    )

    mocked_connect = Mock()
    mocked_mapper = Mock()

    def mocked_reverse_translate(**kwargs):
        return [
            {
                "variant_type": "snv",
                "hgvs_c": f"{kwargs['transcript_accession']}:c.10A>G",
                "hgvs_g": "NC_000001.11:g.100A>G",
            }
        ]

    with (
        patch("src.scripts.reverse_translate_variants.hgvs.dataproviders.uta.connect", return_value=mocked_connect),
        patch("src.scripts.reverse_translate_variants.hgvs.assemblymapper.AssemblyMapper", return_value=mocked_mapper),
        patch("src.scripts.reverse_translate_variants.reverse_translate_hgvs_p", side_effect=mocked_reverse_translate),
    ):
        result = runner.invoke(
            rtv.main,
            [
                "--input",
                str(input_path),
                "--skip",
                "1",
                "--limit",
                "1",
                "--output",
                str(output_path),
            ],
        )

    assert result.exit_code == 0, result.output

    with output_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    assert len(rows) == 1
    assert rows[0]["transcript"] == "NM_B.1"


def test_batch_mode_errors_file_writes_rows_with_error_messages(runner: CliRunner, tmp_path: Path) -> None:
    input_path = tmp_path / "errors_input.tsv"
    output_path = tmp_path / "errors_output.tsv"
    errors_path = tmp_path / "errors.tsv"

    write_tsv(
        input_path,
        rows=[
            {"sample_id": "row1", "transcript": "NM_A.1", "hgvs_p": "p.Arg2His"},
            {"sample_id": "row2", "transcript": "", "hgvs_p": "p.Arg2His"},
            {"sample_id": "row3", "transcript": "NM_B.1", "hgvs_p": "p.Arg2His"},
        ],
        fieldnames=["sample_id", "transcript", "hgvs_p"],
    )

    mocked_connect = Mock()
    mocked_mapper = Mock()

    def mocked_reverse_translate(**kwargs):
        if kwargs["transcript_accession"] == "NM_B.1":
            raise RuntimeError("boom")
        return [
            {
                "variant_type": "snv",
                "hgvs_c": f"{kwargs['transcript_accession']}:c.10A>G",
                "hgvs_g": "NC_000001.11:g.100A>G",
            }
        ]

    with (
        patch("src.scripts.reverse_translate_variants.hgvs.dataproviders.uta.connect", return_value=mocked_connect),
        patch("src.scripts.reverse_translate_variants.hgvs.assemblymapper.AssemblyMapper", return_value=mocked_mapper),
        patch("src.scripts.reverse_translate_variants.reverse_translate_hgvs_p", side_effect=mocked_reverse_translate),
    ):
        result = runner.invoke(
            rtv.main,
            [
                "--input",
                str(input_path),
                "--pass-through-all-columns",
                "--output",
                str(output_path),
                "--errors",
                str(errors_path),
            ],
        )

    assert result.exit_code == 0, result.output

    with output_path.open("r", newline="") as handle:
        output_reader = csv.DictReader(handle, delimiter="\t")
        output_rows = list(output_reader)

    assert len(output_rows) == 3
    assert output_rows[0]["sample_id"] == "row1"
    assert output_rows[0]["variant_type"] == "snv"
    assert output_rows[1]["sample_id"] == "row2"
    assert output_rows[1]["variant_type"] == ""
    assert output_rows[2]["sample_id"] == "row3"
    assert output_rows[2]["variant_type"] == ""

    with errors_path.open("r", newline="") as handle:
        errors_reader = csv.DictReader(handle, delimiter="\t")
        error_rows = list(errors_reader)

    assert errors_reader.fieldnames is not None
    assert "error" in errors_reader.fieldnames
    assert len(error_rows) == 2
    assert error_rows[0]["sample_id"] == "row2"
    assert "missing transcript/UniProt-derived transcript" in error_rows[0]["error"]
    assert error_rows[1]["sample_id"] == "row3"
    assert "failed reverse translation (boom)" in error_rows[1]["error"]
    assert "2 rows written to --errors file" in result.output


def test_single_mode_rejects_skip_option(runner: CliRunner) -> None:
    result = runner.invoke(
        rtv.main,
        [
            "--transcript",
            "NM_TEST.1",
            "--hgvs-p",
            "p.Met1Val",
            "--skip",
            "1",
        ],
    )

    assert result.exit_code != 0
    assert "--skip is only supported with --input." in result.output


def test_single_mode_rejects_errors_option(runner: CliRunner, tmp_path: Path) -> None:
    errors_path = tmp_path / "errors.tsv"

    result = runner.invoke(
        rtv.main,
        [
            "--transcript",
            "NM_TEST.1",
            "--hgvs-p",
            "p.Met1Val",
            "--errors",
            str(errors_path),
        ],
    )

    assert result.exit_code != 0
    assert "--errors is only supported with --input." in result.output


def test_single_mode_one_row_per_input_joins_variants(runner: CliRunner) -> None:
    mocked_connect = Mock()
    mocked_mapper = Mock()

    with (
        patch("src.scripts.reverse_translate_variants.hgvs.dataproviders.uta.connect", return_value=mocked_connect),
        patch("src.scripts.reverse_translate_variants.hgvs.assemblymapper.AssemblyMapper", return_value=mocked_mapper),
        patch(
            "src.scripts.reverse_translate_variants.reverse_translate_hgvs_p",
            return_value=[
                {
                    "variant_type": "snv",
                    "hgvs_c": "NM_TEST.1:c.1A>G",
                    "hgvs_g": "NC_000001.11:g.1A>G",
                },
                {
                    "variant_type": "delins",
                    "hgvs_c": "NM_TEST.1:c.1_3delinsTGG",
                    "hgvs_g": "NC_000001.11:g.1_3delinsTGG",
                },
            ],
        ),
    ):
        result = runner.invoke(
            rtv.main,
            [
                "--transcript",
                "NM_TEST.1",
                "--hgvs-p",
                "p.Ala1Trp",
                "--one-row-per-input",
                "--join-delimiter",
                "|",
            ],
        )

    assert result.exit_code == 0, result.output
    output_lines = [line for line in result.output.splitlines() if line.strip()]
    assert any("hgvs_p\ttranscript\tvariant_type\thgvs_c\thgvs_g" in line for line in output_lines)
    assert any("snv|delins" in line for line in output_lines)
    assert any("NM_TEST.1:c.1A>G|NM_TEST.1:c.1_3delinsTGG" in line for line in output_lines)


def test_batch_mode_one_row_per_input_joins_variants(runner: CliRunner, tmp_path: Path) -> None:
    input_path = tmp_path / "one_row_input.tsv"
    output_path = tmp_path / "one_row_output.tsv"

    write_tsv(
        input_path,
        rows=[
            {"sample_id": "row1", "transcript": "NM_A.1", "hgvs_p": "p.Arg2His"},
            {"sample_id": "row2", "transcript": "NM_B.1", "hgvs_p": "p.Gly2="},
        ],
        fieldnames=["sample_id", "transcript", "hgvs_p"],
    )

    mocked_connect = Mock()
    mocked_mapper = Mock()

    def mocked_reverse_translate(**kwargs):
        if kwargs["transcript_accession"] == "NM_A.1":
            return [
                {
                    "variant_type": "snv",
                    "hgvs_c": "NM_A.1:c.10A>G",
                    "hgvs_g": "NC_000001.11:g.100A>G",
                },
                {
                    "variant_type": "snv",
                    "hgvs_c": "NM_A.1:c.10A>T",
                    "hgvs_g": "NC_000001.11:g.100A>T",
                },
            ]
        return []

    with (
        patch("src.scripts.reverse_translate_variants.hgvs.dataproviders.uta.connect", return_value=mocked_connect),
        patch("src.scripts.reverse_translate_variants.hgvs.assemblymapper.AssemblyMapper", return_value=mocked_mapper),
        patch("src.scripts.reverse_translate_variants.reverse_translate_hgvs_p", side_effect=mocked_reverse_translate),
    ):
        result = runner.invoke(
            rtv.main,
            [
                "--input",
                str(input_path),
                "--pass-through-all-columns",
                "--one-row-per-input",
                "--join-delimiter",
                "|",
                "--output",
                str(output_path),
            ],
        )

    assert result.exit_code == 0, result.output

    with output_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    assert len(rows) == 2
    assert rows[0]["sample_id"] == "row1"
    assert rows[0]["variant_type"] == "snv|snv"
    assert rows[0]["hgvs_c"] == "NM_A.1:c.10A>G|NM_A.1:c.10A>T"
    assert rows[0]["hgvs_g"] == "NC_000001.11:g.100A>G|NC_000001.11:g.100A>T"
    assert len(rows[0]["hgvs_c"].split("|")) == len(rows[0]["hgvs_g"].split("|"))

    assert rows[1]["sample_id"] == "row2"
    assert rows[1]["variant_type"] == ""
    assert rows[1]["hgvs_c"] == ""
    assert rows[1]["hgvs_g"] == ""


def test_single_mode_resolves_uniprot_id(runner: CliRunner) -> None:
    mocked_connect = Mock()
    mocked_mapper = Mock()

    with (
        patch("src.scripts.reverse_translate_variants.hgvs.dataproviders.uta.connect", return_value=mocked_connect),
        patch("src.scripts.reverse_translate_variants.hgvs.assemblymapper.AssemblyMapper", return_value=mocked_mapper),
        patch(
            "src.scripts.reverse_translate_variants.resolve_uniprot_id_to_transcript_accession",
            return_value="NM_TEST.1",
        ) as mocked_resolve_uniprot,
        patch(
            "src.scripts.reverse_translate_variants.reverse_translate_hgvs_p",
            return_value=[
                {
                    "variant_type": "snv",
                    "hgvs_c": "NM_TEST.1:c.1A>G",
                    "hgvs_g": "NC_000001.11:g.1A>G",
                }
            ],
        ) as mocked_reverse_translate,
    ):
        result = runner.invoke(
            rtv.main,
            [
                "--uniprot-id",
                "P04637",
                "--hgvs-p",
                "p.Met1Val",
            ],
        )

    assert result.exit_code == 0, result.output
    assert "NM_TEST.1" in result.output
    assert "NM_TEST.1:c.1A>G" in result.output
    assert mocked_resolve_uniprot.call_count == 1
    assert mocked_reverse_translate.call_count == 1


def test_batch_mode_resolves_unique_uniprot_ids_first(runner: CliRunner, tmp_path: Path) -> None:
    input_path = tmp_path / "input_uniprot.tsv"
    output_path = tmp_path / "output_uniprot.tsv"

    write_tsv(
        input_path,
        rows=[
            {"sample_id": "row1", "uniprot_id": "P04637", "hgvs_p": "p.Arg175His"},
            {"sample_id": "row2", "uniprot_id": "P04637", "hgvs_p": "p.Arg175His"},
            {"sample_id": "row3", "uniprot_id": "P38398", "hgvs_p": "p.Gly2Asp"},
        ],
        fieldnames=["sample_id", "uniprot_id", "hgvs_p"],
    )

    mocked_connect = Mock()
    mocked_mapper = Mock()

    def mocked_resolve_uniprot(**kwargs):
        if kwargs["uniprot_id"] == "P04637":
            return "NM_000546.6"
        if kwargs["uniprot_id"] == "P38398":
            return "NM_000059.4"
        raise AssertionError(f"Unexpected UniProt ID in test: {kwargs['uniprot_id']}")

    def mocked_reverse_translate(**kwargs):
        return [
            {
                "variant_type": "snv",
                "hgvs_c": f"{kwargs['transcript_accession']}:c.1A>G",
                "hgvs_g": "NC_000001.11:g.1A>G",
            }
        ]

    with (
        patch("src.scripts.reverse_translate_variants.hgvs.dataproviders.uta.connect", return_value=mocked_connect),
        patch("src.scripts.reverse_translate_variants.hgvs.assemblymapper.AssemblyMapper", return_value=mocked_mapper),
        patch(
            "src.scripts.reverse_translate_variants.resolve_uniprot_id_to_transcript_accession",
            side_effect=mocked_resolve_uniprot,
        ) as mocked_resolve_uniprot_fn,
        patch("src.scripts.reverse_translate_variants.reverse_translate_hgvs_p", side_effect=mocked_reverse_translate),
    ):
        result = runner.invoke(
            rtv.main,
            [
                "--input",
                str(input_path),
                "--uniprot-column",
                "uniprot_id",
                "--output",
                str(output_path),
            ],
        )

    assert result.exit_code == 0, result.output
    assert mocked_resolve_uniprot_fn.call_count == 2

    with output_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    assert len(rows) == 3
    assert rows[0]["hgvs_c"] == "NM_000546.6:c.1A>G"
    assert rows[1]["hgvs_c"] == "NM_000546.6:c.1A>G"
    assert rows[2]["hgvs_c"] == "NM_000059.4:c.1A>G"

def test_auto_format_hgvs_p_string_with_substitution() -> None:
    """Test auto-formatting single amino acid substitution strings to HGVS p. format."""
    assert rtv.auto_format_hgvs_p_string("A334D") == "p.A334D"
    assert rtv.auto_format_hgvs_p_string("Met1Val") == "p.Met1Val"
    assert rtv.auto_format_hgvs_p_string("*123C") == "p.*123C"
    assert rtv.auto_format_hgvs_p_string("R175=") == "p.R175="


def test_auto_format_hgvs_p_string_with_deletion() -> None:
    """Test auto-formatting single amino acid deletion strings to HGVS p. format."""
    assert rtv.auto_format_hgvs_p_string("A334del") == "p.A334del"
    assert rtv.auto_format_hgvs_p_string("Met1del") == "p.Met1del"
    assert rtv.auto_format_hgvs_p_string("*123del") == "p.*123del"


def test_auto_format_hgvs_p_string_already_formatted() -> None:
    """Test that already formatted HGVS p. strings are returned unchanged."""
    assert rtv.auto_format_hgvs_p_string("p.A334D") == "p.A334D"
    assert rtv.auto_format_hgvs_p_string("p.Met1Val") == "p.Met1Val"
    assert rtv.auto_format_hgvs_p_string("p.A334del") == "p.A334del"


def test_auto_format_hgvs_p_string_non_matching_patterns() -> None:
    """Test that non-matching patterns are returned unchanged."""
    assert rtv.auto_format_hgvs_p_string("invalid") == "invalid"
    assert rtv.auto_format_hgvs_p_string("A334D335E") == "A334D335E"  # Multiple changes
    assert rtv.auto_format_hgvs_p_string("p.Arg175His") == "p.Arg175His"  # Already proper format
    assert rtv.auto_format_hgvs_p_string("") == ""  # Empty string


def test_auto_format_hgvs_p_string_with_whitespace() -> None:
    """Test that whitespace is trimmed before formatting."""
    assert rtv.auto_format_hgvs_p_string("  A334D  ") == "p.A334D"
    assert rtv.auto_format_hgvs_p_string("  A334del  ") == "p.A334del"


def test_single_mode_auto_format_hgvs_p_substitution(runner: CliRunner, tmp_path: Path) -> None:
    """Test --auto-format-hgvs-p flag in single mode with amino acid substitution."""
    output_path = tmp_path / "output.tsv"

    data_provider = Mock()
    data_provider.get_tx_identity_info.return_value = {"cds_start_i": 0, "cds_end_i": 3}
    data_provider.get_seq.return_value = "ATG"  # Single codon = M (Met)

    mocked_mapper = Mock()

    with patch(
        "src.scripts.reverse_translate_variants.hgvs.dataproviders.uta.connect",
        return_value=data_provider,
    ), patch(
        "src.scripts.reverse_translate_variants.hgvs.assemblymapper.AssemblyMapper",
        return_value=mocked_mapper,
    ), patch(
        "src.scripts.reverse_translate_variants.map_hgvs_c_to_hgvs_g",
        side_effect=lambda parser, mapper, tx, hgvs_c: f"GENOMIC:{tx}:{hgvs_c}",
    ):
        result = runner.invoke(
            rtv.main,
            [
                "--transcript",
                "NM_TEST.1",
                "--hgvs-p",
                "M1V",  # One-letter format without p. prefix, position 1 is valid for 1 AA sequence
                "--auto-format-hgvs-p",
                "--output",
                str(output_path),
            ],
        )

    assert result.exit_code == 0, result.output
    assert "Generated" in result.output


def test_single_mode_auto_format_hgvs_p_deletion(runner: CliRunner, tmp_path: Path) -> None:
    """Test --auto-format-hgvs-p flag in single mode with amino acid deletion."""
    output_path = tmp_path / "output.tsv"

    data_provider = Mock()
    data_provider.get_tx_identity_info.return_value = {"cds_start_i": 0, "cds_end_i": 6}
    data_provider.get_seq.return_value = "ATGGCT"  # Two codons = M (Met), A (Ala)

    mocked_mapper = Mock()

    with patch(
        "src.scripts.reverse_translate_variants.hgvs.dataproviders.uta.connect",
        return_value=data_provider,
    ), patch(
        "src.scripts.reverse_translate_variants.hgvs.assemblymapper.AssemblyMapper",
        return_value=mocked_mapper,
    ), patch(
        "src.scripts.reverse_translate_variants.map_hgvs_c_to_hgvs_g",
        side_effect=lambda parser, mapper, tx, hgvs_c: f"GENOMIC:{tx}:{hgvs_c}",
    ), patch(
        "src.scripts.reverse_translate_variants.reverse_translate_hgvs_p",
        return_value=[  # Mock a successful deletion result
            {
                "variant_type": "del",
                "hgvs_c": "NM_TEST.1:c.1_3del",
                "hgvs_g": "NC_000001.11:g.1_3del",
            }
        ],
    ):
        result = runner.invoke(
            rtv.main,
            [
                "--transcript",
                "NM_TEST.1",
                "--hgvs-p",
                "M1del",  # Deletion format without p. prefix
                "--auto-format-hgvs-p",
                "--output",
                str(output_path),
            ],
        )

    assert result.exit_code == 0, result.output
    assert "Generated" in result.output


def test_batch_mode_auto_format_hgvs_p(runner: CliRunner, tmp_path: Path) -> None:
    """Test --auto-format-hgvs-p flag in batch mode with mixed substitutions and deletions."""
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"

    write_tsv(
        input_path,
        rows=[
            {"sample_id": "row1", "transcript": "NM_A.1", "hgvs_p": "A334D"},  # Substitution without p.
            {"sample_id": "row2", "transcript": "NM_B.1", "hgvs_p": "M1V"},  # Short form
            {"sample_id": "row3", "transcript": "NM_C.1", "hgvs_p": "A334del"},  # Deletion without p.
        ],
        fieldnames=["sample_id", "transcript", "hgvs_p"],
    )

    mocked_connect = Mock()
    mocked_mapper = Mock()

    def mocked_reverse_translate(**kwargs):
        return [
            {
                "variant_type": "snv",
                "hgvs_c": f"{kwargs['transcript_accession']}:c.1A>G",
                "hgvs_g": "NC_000001.11:g.1A>G",
            }
        ]

    with (
        patch("src.scripts.reverse_translate_variants.hgvs.dataproviders.uta.connect", return_value=mocked_connect),
        patch("src.scripts.reverse_translate_variants.hgvs.assemblymapper.AssemblyMapper", return_value=mocked_mapper),
        patch("src.scripts.reverse_translate_variants.reverse_translate_hgvs_p", side_effect=mocked_reverse_translate),
    ):
        result = runner.invoke(
            rtv.main,
            [
                "--input",
                str(input_path),
                "--auto-format-hgvs-p",
                "--output",
                str(output_path),
            ],
        )

    assert result.exit_code == 0, result.output

    with output_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    assert len(rows) == 3
    assert rows[0]["hgvs_c"] == "NM_A.1:c.1A>G"
    assert rows[1]["hgvs_c"] == "NM_B.1:c.1A>G"
    assert rows[2]["hgvs_c"] == "NM_C.1:c.1A>G"


def test_parse_hgvs_protein_deletion() -> None:
    """Test parsing HGVS p. strings with deletion suffix."""
    parsed = rtv.parse_hgvs_protein_change("p.Arg175del")
    assert parsed.reference_aa == "R"
    assert parsed.position == 175
    assert parsed.alternate_aa == ""

    parsed_single_letter = rtv.parse_hgvs_protein_change("p.A334del")
    assert parsed_single_letter.reference_aa == "A"
    assert parsed_single_letter.position == 334
    assert parsed_single_letter.alternate_aa == ""

    parsed_with_accession = rtv.parse_hgvs_protein_change("NP_000537.3:p.Met1del")
    assert parsed_with_accession.reference_aa == "M"
    assert parsed_with_accession.position == 1
    assert parsed_with_accession.alternate_aa == ""


def test_auto_format_hgvs_p_string_with_dash_deletion() -> None:
    """Test auto-formatting deletion strings with dash notation."""
    assert rtv.auto_format_hgvs_p_string("A334-") == "p.A334del"
    assert rtv.auto_format_hgvs_p_string("M1-") == "p.M1del"
    assert rtv.auto_format_hgvs_p_string("  A334-  ") == "p.A334del"


def test_single_mode_deletion_variant(runner: CliRunner, tmp_path: Path) -> None:
    """Test reverse translation of a single amino acid deletion."""
    output_path = tmp_path / "output.tsv"

    data_provider = Mock()
    data_provider.get_tx_identity_info.return_value = {"cds_start_i": 0, "cds_end_i": 9}
    data_provider.get_seq.return_value = "ATGGCTAAA"  # Three codons = M, A, K

    mocked_mapper = Mock()

    with patch(
        "src.scripts.reverse_translate_variants.hgvs.dataproviders.uta.connect",
        return_value=data_provider,
    ), patch(
        "src.scripts.reverse_translate_variants.hgvs.assemblymapper.AssemblyMapper",
        return_value=mocked_mapper,
    ), patch(
        "src.scripts.reverse_translate_variants.map_hgvs_c_to_hgvs_g",
        side_effect=lambda parser, mapper, tx, hgvs_c: f"GENOMIC:{tx}:{hgvs_c}",
    ):
        result = runner.invoke(
            rtv.main,
            [
                "--transcript",
                "NM_TEST.1",
                "--hgvs-p",
                "p.A2del",  # Delete the Ala at position 2
                "--output",
                str(output_path),
            ],
        )

    assert result.exit_code == 0, result.output
    
    # Verify output contains a deletion variant
    with output_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    assert len(rows) == 1
    assert rows[0]["variant_type"] == "del"
    assert "del" in rows[0]["hgvs_c"]


def test_batch_mode_deletion_with_auto_format(runner: CliRunner, tmp_path: Path) -> None:
    """Test batch mode with deletion auto-format from dash notation."""
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"

    write_tsv(
        input_path,
        rows=[
            {"id": "var1", "transcript": "NM_A.1", "hgvs_p": "M1-"},  # Deletion with dash
            {"id": "var2", "transcript": "NM_B.1", "hgvs_p": "A2del"},  # Deletion already proper format
        ],
        fieldnames=["id", "transcript", "hgvs_p"],
    )

    mocked_connect = Mock()
    mocked_mapper = Mock()

    def mocked_reverse_translate(**kwargs):
        return [
            {
                "variant_type": "del",
                "hgvs_c": f"{kwargs['transcript_accession']}:c.1_3del",
                "hgvs_g": "NC_000001.11:g.1_3del",
            }
        ]

    with (
        patch("src.scripts.reverse_translate_variants.hgvs.dataproviders.uta.connect", return_value=mocked_connect),
        patch("src.scripts.reverse_translate_variants.hgvs.assemblymapper.AssemblyMapper", return_value=mocked_mapper),
        patch("src.scripts.reverse_translate_variants.reverse_translate_hgvs_p", side_effect=mocked_reverse_translate),
    ):
        result = runner.invoke(
            rtv.main,
            [
                "--input",
                str(input_path),
                "--auto-format-hgvs-p",
                "--output",
                str(output_path),
            ],
        )

    assert result.exit_code == 0, result.output

    with output_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    assert len(rows) == 2
    assert rows[0]["variant_type"] == "del"
    assert rows[1]["variant_type"] == "del"


def test_batch_mode_accepts_csv_input_format(runner: CliRunner, tmp_path: Path) -> None:
    input_path = tmp_path / "input.csv"
    output_path = tmp_path / "output.tsv"

    write_csv(
        input_path,
        rows=[
            {"sample_id": "row1", "transcript": "NM_A.1", "hgvs_p": "p.Arg2His"},
            {"sample_id": "row2", "transcript": "NM_B.1", "hgvs_p": "p.Gly2="},
        ],
        fieldnames=["sample_id", "transcript", "hgvs_p"],
    )

    mocked_connect = Mock()
    mocked_mapper = Mock()

    def mocked_reverse_translate(**kwargs):
        if kwargs["transcript_accession"] == "NM_A.1":
            return [
                {
                    "variant_type": "snv",
                    "hgvs_c": "NM_A.1:c.10A>G",
                    "hgvs_g": "NC_000001.11:g.100A>G",
                }
            ]
        return []

    with (
        patch("src.scripts.reverse_translate_variants.hgvs.dataproviders.uta.connect", return_value=mocked_connect),
        patch("src.scripts.reverse_translate_variants.hgvs.assemblymapper.AssemblyMapper", return_value=mocked_mapper),
        patch("src.scripts.reverse_translate_variants.reverse_translate_hgvs_p", side_effect=mocked_reverse_translate),
    ):
        result = runner.invoke(
            rtv.main,
            [
                "--input",
                str(input_path),
                "--input-format",
                "csv",
                "--output",
                str(output_path),
            ],
        )

    assert result.exit_code == 0, result.output

    with output_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    assert len(rows) == 2
    assert rows[0]["hgvs_c"] == "NM_A.1:c.10A>G"
    assert rows[0]["transcript"] == "NM_A.1"
    assert rows[1]["transcript"] == "NM_B.1"
    assert rows[1]["variant_type"] == ""


def test_batch_mode_honors_csv_field_size_limit_option(runner: CliRunner, tmp_path: Path) -> None:
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"

    write_tsv(
        input_path,
        rows=[
            {
                "transcript": "NM_A.1",
                "hgvs_p": "p.Arg2His",
                "note": "X" * 200,
            }
        ],
        fieldnames=["transcript", "hgvs_p", "note"],
    )

    with (
        patch("src.scripts.reverse_translate_variants.hgvs.dataproviders.uta.connect", return_value=Mock()),
        patch("src.scripts.reverse_translate_variants.hgvs.assemblymapper.AssemblyMapper", return_value=Mock()),
        patch(
            "src.scripts.reverse_translate_variants.reverse_translate_hgvs_p",
            return_value=[
                {
                    "variant_type": "snv",
                    "hgvs_c": "NM_A.1:c.10A>G",
                    "hgvs_g": "NC_000001.11:g.100A>G",
                }
            ],
        ),
    ):
        result = runner.invoke(
            rtv.main,
            [
                "--input",
                str(input_path),
                "--csv-field-size-limit",
                "50",
                "--output",
                str(output_path),
            ],
        )

    assert result.exit_code != 0
    assert isinstance(result.exception, csv.Error)
    assert "field larger than field limit" in str(result.exception)
