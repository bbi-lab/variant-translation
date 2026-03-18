"""
Reverse-translate a protein HGVS change to DNA HGVS variants.

Given a transcript accession and a protein HGVS p. change (for example, `p.Gly12Asp`), this script enumerates all
single-nucleotide variants in the affected codon that yield the same amino-acid change. Optionally, it can also
consider small codon-local indels (up to a configurable size, default 3 nt) and keep only those that reproduce the
same protein-level change.

Output is tab-separated with both HGVS c. and HGVS g. expressions.

Examples:
  python -m src.scripts.reverse_translate_variants \
	--transcript NM_000546.6 \
	--hgvs-p p.Arg175His

  python -m src.scripts.reverse_translate_variants \
	--uniprot-id P04637 \
	--hgvs-p p.Arg175His

  python -m src.scripts.reverse_translate_variants \
	--transcript NM_000546.6 \
	--hgvs-p p.Arg175His \
	--include-indels \
	--output tp53_r175h_reverse_translation.tsv

  python -m src.scripts.reverse_translate_variants \
	--input input.tsv \
	--uniprot-column uniprot_id \
	--hgvs-p-column hgvs_p \
	--include-indels \
	--one-row-per-input \
	--join-delimiter "|" \
	--output reverse_translation_joined.tsv

  python -m src.scripts.reverse_translate_variants \
	--input input.tsv \
	--pass-through-all-columns \
	--pass-through-prefix input_ \
	--output reverse_translation_with_input_columns.tsv
"""

import csv
from dataclasses import dataclass
from itertools import product
import json
from pathlib import Path
import re
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

import click
from dotenv import load_dotenv
import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.parser

load_dotenv()

DNA_BASES = ("A", "C", "G", "T")
DNA_COMPLEMENTS = str.maketrans({"A": "T", "C": "G", "G": "C", "T": "A"})


CODON_TO_AA = {
	"TTT": "F",
	"TTC": "F",
	"TTA": "L",
	"TTG": "L",
	"TCT": "S",
	"TCC": "S",
	"TCA": "S",
	"TCG": "S",
	"TAT": "Y",
	"TAC": "Y",
	"TAA": "*",
	"TAG": "*",
	"TGT": "C",
	"TGC": "C",
	"TGA": "*",
	"TGG": "W",
	"CTT": "L",
	"CTC": "L",
	"CTA": "L",
	"CTG": "L",
	"CCT": "P",
	"CCC": "P",
	"CCA": "P",
	"CCG": "P",
	"CAT": "H",
	"CAC": "H",
	"CAA": "Q",
	"CAG": "Q",
	"CGT": "R",
	"CGC": "R",
	"CGA": "R",
	"CGG": "R",
	"ATT": "I",
	"ATC": "I",
	"ATA": "I",
	"ATG": "M",
	"ACT": "T",
	"ACC": "T",
	"ACA": "T",
	"ACG": "T",
	"AAT": "N",
	"AAC": "N",
	"AAA": "K",
	"AAG": "K",
	"AGT": "S",
	"AGC": "S",
	"AGA": "R",
	"AGG": "R",
	"GTT": "V",
	"GTC": "V",
	"GTA": "V",
	"GTG": "V",
	"GCT": "A",
	"GCC": "A",
	"GCA": "A",
	"GCG": "A",
	"GAT": "D",
	"GAC": "D",
	"GAA": "E",
	"GAG": "E",
	"GGT": "G",
	"GGC": "G",
	"GGA": "G",
	"GGG": "G",
}


AA_3_TO_1 = {
	"Ala": "A",
	"Arg": "R",
	"Asn": "N",
	"Asp": "D",
	"Cys": "C",
	"Gln": "Q",
	"Glu": "E",
	"Gly": "G",
	"His": "H",
	"Ile": "I",
	"Leu": "L",
	"Lys": "K",
	"Met": "M",
	"Phe": "F",
	"Pro": "P",
	"Ser": "S",
	"Thr": "T",
	"Trp": "W",
	"Tyr": "Y",
	"Val": "V",
	"Ter": "*",
}


HGVS_PROTEIN_RE = re.compile(
	r"^(?:(?P<ac>[^:]+):)?p\.\(?(?P<ref>[A-Za-z\*]{1,4})(?P<pos>\d+)(?P<alt>[A-Za-z\*=]{1,4})\)?$"
)


UNIPROT_MANE_API_TEMPLATE = "https://rest.uniprot.org/uniprotkb/{uniprot_id}?fields=xref_mane-select"


@dataclass(frozen=True)
class ProteinChange:
	reference_aa: str
	position: int
	alternate_aa: str


@dataclass(frozen=True)
class CandidateVariant:
	variant_type: str
	hgvs_c: str


@dataclass(frozen=True)
class ManeSelectCrossReference:
	ensembl_transcript_id: str
	refseq_protein_id: str | None
	refseq_nucleotide_id: str | None


def aa_token_to_one_letter(token: str) -> str:
	normalized = token.strip()
	if normalized == "*":
		return "*"

	upper = normalized.upper()
	if upper in {"TER", "STOP"}:
		return "*"

	if len(normalized) == 1 and normalized.isalpha():
		return upper

	if len(normalized) == 3:
		title = normalized[0].upper() + normalized[1:].lower()
		if title in AA_3_TO_1:
			return AA_3_TO_1[title]

	raise click.ClickException(f"Unsupported amino acid token: {token}")


def parse_hgvs_protein_change(hgvs_p: str) -> ProteinChange:
	hgvs_normalized = hgvs_p.strip()
	
	deletion_pattern = re.compile(r"^(?:(?P<ac>[^:]+):)?p\.\(?(?P<ref>[A-Za-z\*]{1,4})(?P<pos>\d+)del\)?$")
	deletion_match = deletion_pattern.match(hgvs_normalized)
	if deletion_match is not None:
		ref_token = deletion_match.group("ref")
		position = int(deletion_match.group("pos"))
		reference_aa = aa_token_to_one_letter(ref_token)
		return ProteinChange(reference_aa=reference_aa, position=position, alternate_aa="")
	
	match = HGVS_PROTEIN_RE.match(hgvs_normalized)
	if match is None:
		raise click.ClickException(
			"Could not parse HGVS protein string. Expected forms like 'p.Arg175His', 'p.Arg175del', or 'NP_000537.3:p.Arg175His'."
		)

	ref_token = match.group("ref")
	alt_token = match.group("alt")
	position = int(match.group("pos"))

	reference_aa = aa_token_to_one_letter(ref_token)
	alternate_aa = reference_aa if alt_token == "=" else aa_token_to_one_letter(alt_token)

	return ProteinChange(reference_aa=reference_aa, position=position, alternate_aa=alternate_aa)


def translate_cds(coding_sequence: str) -> str:
	amino_acids: list[str] = []
	coding_length = len(coding_sequence) - (len(coding_sequence) % 3)
	for codon_start in range(0, coding_length, 3):
		codon = coding_sequence[codon_start : codon_start + 3]
		amino_acid = CODON_TO_AA.get(codon, "X")
		amino_acids.append(amino_acid)
		if amino_acid == "*":
			break
	return "".join(amino_acids)


def apply_cds_edit(coding_sequence: str, edit_start_index: int, deleted_length: int, inserted_sequence: str) -> str:
	return (
		coding_sequence[:edit_start_index]
		+ inserted_sequence
		+ coding_sequence[edit_start_index + deleted_length :]
	)


def matches_requested_protein_change(
	reference_protein: str,
	alternate_protein: str,
	protein_change: ProteinChange,
) -> bool:
	aa_position = protein_change.position
	
	if protein_change.alternate_aa == "":
		if aa_position > len(reference_protein):
			return False
		position_index = aa_position - 1
		if reference_protein[position_index] != protein_change.reference_aa:
			return False
		if len(alternate_protein) + 1 != len(reference_protein):
			return False
		if reference_protein[:position_index] != alternate_protein[:position_index]:
			return False
		if reference_protein[position_index + 1 :] != alternate_protein[position_index:]:
			return False
		return True
	
	if aa_position > len(reference_protein) or aa_position > len(alternate_protein):
		return False

	position_index = aa_position - 1

	if reference_protein[position_index] != protein_change.reference_aa:
		return False

	if alternate_protein[position_index] != protein_change.alternate_aa:
		return False

	if reference_protein[:position_index] != alternate_protein[:position_index]:
		return False

	if protein_change.alternate_aa == "*":
		return True

	return reference_protein[position_index + 1 :] == alternate_protein[position_index + 1 :]


def hgvs_c_for_substitution(codon_c_start: int, codon_offset: int, reference_nt: str, alternate_nt: str) -> str:
	c_position = codon_c_start + codon_offset
	return f"c.{c_position}{reference_nt}>{alternate_nt}"


def hgvs_c_for_indel(codon_c_start: int, codon_offset: int, deleted_length: int, inserted_sequence: str) -> str:
	if deleted_length == 0:
		left_position = codon_c_start + codon_offset - 1
		right_position = left_position + 1
		return f"c.{left_position}_{right_position}ins{inserted_sequence}"

	start_position = codon_c_start + codon_offset
	end_position = start_position + deleted_length - 1

	if not inserted_sequence:
		if deleted_length == 1:
			return f"c.{start_position}del"
		return f"c.{start_position}_{end_position}del"

	if deleted_length == 1:
		return f"c.{start_position}delins{inserted_sequence}"
	return f"c.{start_position}_{end_position}delins{inserted_sequence}"


def reverse_complement_dna(sequence: str) -> str:
	return sequence.upper().translate(DNA_COMPLEMENTS)[::-1]


def is_inversion_replacement(deleted_sequence: str, inserted_sequence: str) -> bool:
	return (
		len(deleted_sequence) > 1
		and len(deleted_sequence) == len(inserted_sequence)
		and inserted_sequence == reverse_complement_dna(deleted_sequence)
	)


def minimize_same_length_delins(
	codon_offset: int,
	deleted_sequence: str,
	inserted_sequence: str,
) -> tuple[int, str, str] | None:
	if len(deleted_sequence) != len(inserted_sequence):
		return codon_offset, deleted_sequence, inserted_sequence

	leading_matches = 0
	while (
		leading_matches < len(deleted_sequence)
		and deleted_sequence[leading_matches] == inserted_sequence[leading_matches]
	):
		leading_matches += 1

	if leading_matches == len(deleted_sequence):
		return None

	trailing_matches = 0
	while (
		trailing_matches < len(deleted_sequence) - leading_matches
		and deleted_sequence[-(trailing_matches + 1)] == inserted_sequence[-(trailing_matches + 1)]
	):
		trailing_matches += 1

	trimmed_deleted_sequence = deleted_sequence[leading_matches : len(deleted_sequence) - trailing_matches]
	trimmed_inserted_sequence = inserted_sequence[leading_matches : len(inserted_sequence) - trailing_matches]

	if len(trimmed_deleted_sequence) == 1:
		return None

	return codon_offset + leading_matches, trimmed_deleted_sequence, trimmed_inserted_sequence


def is_frame_preserving_indel(deleted_length: int, inserted_length: int) -> bool:
	return (inserted_length - deleted_length) % 3 == 0


def hgvs_c_for_delins_or_inv(
	codon_c_start: int,
	codon_offset: int,
	deleted_sequence: str,
	inserted_sequence: str,
	use_inv_notation: bool,
) -> tuple[str, str]:
	deleted_length = len(deleted_sequence)
	if use_inv_notation and is_inversion_replacement(deleted_sequence, inserted_sequence):
		start_position = codon_c_start + codon_offset
		end_position = start_position + deleted_length - 1
		return "inv", f"c.{start_position}_{end_position}inv"

	return "delins", hgvs_c_for_indel(codon_c_start, codon_offset, deleted_length, inserted_sequence)


def enumerate_snv_candidates(
	coding_sequence: str,
	codon_sequence: str,
	codon_c_start: int,
	codon_cds_start_index: int,
	reference_protein: str,
	protein_change: ProteinChange,
) -> list[CandidateVariant]:
	candidates: list[CandidateVariant] = []

	for codon_offset, reference_nt in enumerate(codon_sequence):
		for alternate_nt in DNA_BASES:
			if alternate_nt == reference_nt:
				continue

			alternate_cds = apply_cds_edit(
				coding_sequence,
				codon_cds_start_index + codon_offset,
				1,
				alternate_nt,
			)
			alternate_protein = translate_cds(alternate_cds)
			if not matches_requested_protein_change(reference_protein, alternate_protein, protein_change):
				continue

			hgvs_c = hgvs_c_for_substitution(codon_c_start, codon_offset, reference_nt, alternate_nt)
			candidates.append(CandidateVariant(variant_type="snv", hgvs_c=hgvs_c))

	return candidates


def enumerate_indel_candidates(
	coding_sequence: str,
	codon_sequence: str,
	codon_c_start: int,
	codon_cds_start_index: int,
	reference_protein: str,
	protein_change: ProteinChange,
	max_indel_size: int,
	use_inv_notation: bool,
) -> list[CandidateVariant]:
	candidates: list[CandidateVariant] = []

	# Pure insertions (internal to codon only): c.N_N+1insX
	for codon_offset in (1, 2):
		for inserted_length in range(1, max_indel_size + 1):
			if not is_frame_preserving_indel(deleted_length=0, inserted_length=inserted_length):
				continue
			for inserted_bases in product(DNA_BASES, repeat=inserted_length):
				inserted_sequence = "".join(inserted_bases)
				alternate_cds = apply_cds_edit(
					coding_sequence,
					codon_cds_start_index + codon_offset,
					0,
					inserted_sequence,
				)
				alternate_protein = translate_cds(alternate_cds)
				if not matches_requested_protein_change(reference_protein, alternate_protein, protein_change):
					continue

				hgvs_c = hgvs_c_for_indel(codon_c_start, codon_offset, 0, inserted_sequence)
				candidates.append(CandidateVariant(variant_type="insertion", hgvs_c=hgvs_c))

	# Deletions and delins fully contained in the codon.
	for deleted_length in range(1, min(3, max_indel_size) + 1):
		for codon_offset in range(0, 4 - deleted_length):
			deleted_sequence = coding_sequence[
				codon_cds_start_index + codon_offset : codon_cds_start_index + codon_offset + deleted_length
			]

			# Pure deletion
			if is_frame_preserving_indel(deleted_length=deleted_length, inserted_length=0):
				alternate_cds = apply_cds_edit(
					coding_sequence,
					codon_cds_start_index + codon_offset,
					deleted_length,
					"",
				)
				alternate_protein = translate_cds(alternate_cds)
				if matches_requested_protein_change(reference_protein, alternate_protein, protein_change):
					hgvs_c = hgvs_c_for_indel(codon_c_start, codon_offset, deleted_length, "")
					candidates.append(CandidateVariant(variant_type="deletion", hgvs_c=hgvs_c))

			# Delins
			for inserted_length in range(1, max_indel_size + 1):
				if not is_frame_preserving_indel(deleted_length=deleted_length, inserted_length=inserted_length):
					continue
				for inserted_bases in product(DNA_BASES, repeat=inserted_length):
					inserted_sequence = "".join(inserted_bases)
					trimmed_indel = minimize_same_length_delins(
						codon_offset,
						deleted_sequence,
						inserted_sequence,
					)
					if inserted_length == deleted_length:
						if trimmed_indel is None:
							continue
						trimmed_codon_offset, trimmed_deleted_sequence, trimmed_inserted_sequence = trimmed_indel
					else:
						trimmed_codon_offset = codon_offset
						trimmed_deleted_sequence = deleted_sequence
						trimmed_inserted_sequence = inserted_sequence

					alternate_cds = apply_cds_edit(
						coding_sequence,
						codon_cds_start_index + codon_offset,
						deleted_length,
						inserted_sequence,
					)
					alternate_protein = translate_cds(alternate_cds)
					if not matches_requested_protein_change(reference_protein, alternate_protein, protein_change):
						continue

					variant_type, hgvs_c = hgvs_c_for_delins_or_inv(
						codon_c_start,
						trimmed_codon_offset,
						trimmed_deleted_sequence,
						trimmed_inserted_sequence,
						use_inv_notation=use_inv_notation,
					)
					candidates.append(CandidateVariant(variant_type=variant_type, hgvs_c=hgvs_c))

	return candidates


def map_hgvs_c_to_hgvs_g(
	parser: hgvs.parser.Parser,
	mapper: hgvs.assemblymapper.AssemblyMapper,
	transcript_accession: str,
	hgvs_c: str,
) -> str:
	variant_c = parser.parse_hgvs_variant(f"{transcript_accession}:{hgvs_c}")
	variant_g = mapper.c_to_g(variant_c)
	return str(variant_g)


def looks_like_transcript_accession(accession: str) -> bool:
	return bool(re.match(r"^(?:NM_|NR_|XM_|XR_|ENST|LRG_)", accession))


def iter_string_values(value: Any) -> Any:
	if isinstance(value, str):
		yield value
		return

	if isinstance(value, dict):
		for nested_value in value.values():
			yield from iter_string_values(nested_value)
		return

	if isinstance(value, (list, tuple)):
		for nested_value in value:
			yield from iter_string_values(nested_value)


def first_transcript_accession_from_protein_mapping_result(mapping_result: Any) -> str | None:
	if mapping_result is None:
		return None

	for candidate in iter_string_values(mapping_result):
		if looks_like_transcript_accession(candidate):
			return candidate

	return None


def resolve_transcript_from_refseq_protein_id(data_provider: Any, refseq_protein_id: str) -> str | None:
	if not hasattr(data_provider, "get_tx_for_pro_ac"):
		return None

	try:
		mapping_result = data_provider.get_tx_for_pro_ac(refseq_protein_id)
	except Exception:
		return None

	return first_transcript_accession_from_protein_mapping_result(mapping_result)


def fetch_uniprot_mane_select_cross_reference(uniprot_id: str) -> ManeSelectCrossReference:
	request_url = UNIPROT_MANE_API_TEMPLATE.format(uniprot_id=uniprot_id)

	try:
		with urlopen(request_url, timeout=30) as response:
			response_payload = json.loads(response.read().decode("utf-8"))
	except HTTPError as exception:
		raise click.ClickException(f"UniProt lookup failed for {uniprot_id}: HTTP {exception.code}") from exception
	except URLError as exception:
		raise click.ClickException(f"UniProt lookup failed for {uniprot_id}: {exception.reason}") from exception
	except json.JSONDecodeError as exception:
		raise click.ClickException(f"UniProt lookup returned invalid JSON for {uniprot_id}.") from exception

	cross_references = response_payload.get("uniProtKBCrossReferences") or []
	if not cross_references:
		raise click.ClickException(f"No MANE Select cross-reference found in UniProt for {uniprot_id}.")

	first_cross_reference = cross_references[0]
	ensembl_transcript_id = (first_cross_reference.get("id") or "").strip()
	if not ensembl_transcript_id:
		raise click.ClickException(f"UniProt MANE Select data for {uniprot_id} is missing cross-reference id.")

	refseq_protein_id: str | None = None
	refseq_nucleotide_id: str | None = None

	for property_item in first_cross_reference.get("properties") or []:
		property_key = (property_item.get("key") or "").strip()
		property_value = (property_item.get("value") or "").strip()
		if not property_value:
			continue

		if property_key == "RefSeqProteinId":
			refseq_protein_id = property_value
		elif property_key in {"RefSeqNucleotideId", "RefSeqNucleotided"}:
			refseq_nucleotide_id = property_value

	return ManeSelectCrossReference(
		ensembl_transcript_id=ensembl_transcript_id,
		refseq_protein_id=refseq_protein_id,
		refseq_nucleotide_id=refseq_nucleotide_id,
	)


def resolve_uniprot_id_to_transcript_accession(
	uniprot_id: str,
	uniprot_target: str,
	data_provider: Any,
	uniprot_cache: dict[str, str],
) -> str:
	normalized_uniprot_id = uniprot_id.strip()
	if not normalized_uniprot_id:
		raise click.ClickException("UniProt ID cannot be empty.")

	cache_key = f"{uniprot_target}:{normalized_uniprot_id}"
	if cache_key in uniprot_cache:
		return uniprot_cache[cache_key]

	cross_reference = fetch_uniprot_mane_select_cross_reference(normalized_uniprot_id)

	if uniprot_target == "ensembl-transcript":
		transcript_accession = cross_reference.ensembl_transcript_id
	elif uniprot_target == "refseq-protein":
		if cross_reference.refseq_protein_id is None:
			raise click.ClickException(
				f"UniProt MANE Select data for {normalized_uniprot_id} does not include RefSeqProteinId."
			)

		transcript_accession = resolve_transcript_from_refseq_protein_id(
			data_provider=data_provider,
			refseq_protein_id=cross_reference.refseq_protein_id,
		) or (cross_reference.refseq_nucleotide_id or "")

		if not transcript_accession:
			raise click.ClickException(
				f"Unable to resolve a transcript accession for UniProt {normalized_uniprot_id} "
				f"from RefSeq protein {cross_reference.refseq_protein_id}."
			)
	else:
		raise click.ClickException(f"Unsupported --uniprot-target value: {uniprot_target}")

	uniprot_cache[cache_key] = transcript_accession
	return transcript_accession


def auto_format_hgvs_p_string(hgvs_p: str) -> str:
	"""
	Auto-format protein HGVS strings from one-letter notation (e.g., 'A334D' or 'A334del') to HGVS p. format.

	If the string already starts with 'p.', it is returned as-is.
	If it matches a single-amino-acid substitution pattern (e.g., 'A334D', 'Met1Val', 'R175=') or deletion
	pattern (e.g., 'A334del', 'Met1del', 'A334-'), the 'p.' prefix is prepended.
	Otherwise, the string is returned unchanged.

	Args:
		hgvs_p: Input string that may or may not be in HGVS p. format.

	Returns:
		Formatted string with 'p.' prefix if it matches single-AA patterns.
	"""
	hgvs_p = hgvs_p.strip()

	if not hgvs_p:
		return hgvs_p

	if hgvs_p.startswith("p."):
		return hgvs_p

	single_aa_substitution_pattern = r"^[A-Z\*][a-z]{0,2}\d+([A-Z\*][a-z]{0,2}|=)$"
	single_aa_deletion_pattern = r"^[A-Z\*][a-z]{0,2}\d+(?:del|-)$"

	if re.match(single_aa_substitution_pattern, hgvs_p) or re.match(single_aa_deletion_pattern, hgvs_p):
		formatted = f"p.{hgvs_p}"
		if formatted.endswith("-"):
			formatted = formatted[:-1] + "del"
		return formatted

	return hgvs_p


def get_coding_sequence_and_reference_protein(
	data_provider: Any,
	transcript_accession: str,
	transcript_cache: dict[str, tuple[str, str]],
) -> tuple[str, str]:
	if transcript_accession in transcript_cache:
		return transcript_cache[transcript_accession]

	transcript_info = data_provider.get_tx_identity_info(transcript_accession)
	if transcript_info is None:
		raise click.ClickException(f"Unable to fetch transcript metadata for: {transcript_accession}")

	full_transcript_sequence = data_provider.get_seq(transcript_accession)
	if not full_transcript_sequence:
		raise click.ClickException(f"Unable to fetch transcript sequence for: {transcript_accession}")

	cds_start_index = int(transcript_info["cds_start_i"])
	cds_end_index = int(transcript_info["cds_end_i"])
	coding_sequence = full_transcript_sequence[cds_start_index:cds_end_index].upper()
	reference_protein = translate_cds(coding_sequence)

	transcript_cache[transcript_accession] = (coding_sequence, reference_protein)
	return coding_sequence, reference_protein


def reverse_translate_hgvs_p(
	transcript_accession: str,
	hgvs_protein: str,
	include_indels: bool,
	max_indel_size: int,
	strict_ref_aa: bool,
	use_inv_notation: bool,
	parser: hgvs.parser.Parser,
	mapper: hgvs.assemblymapper.AssemblyMapper,
	data_provider: Any,
	transcript_cache: dict[str, tuple[str, str]],
) -> list[dict[str, str]]:
	protein_change = parse_hgvs_protein_change(hgvs_protein)

	if protein_change.reference_aa == "*":
		raise click.ClickException("Stop-loss reverse translation is not supported by this script.")
	
	is_deletion = protein_change.alternate_aa == ""

	coding_sequence, reference_protein = get_coding_sequence_and_reference_protein(
		data_provider,
		transcript_accession,
		transcript_cache,
	)

	codon_cds_start_index = (protein_change.position - 1) * 3
	codon_cds_end_index = codon_cds_start_index + 3
	if codon_cds_end_index > len(coding_sequence):
		raise click.ClickException(
			f"Protein position {protein_change.position} exceeds coding sequence length for {transcript_accession}."
		)

	codon_sequence = coding_sequence[codon_cds_start_index:codon_cds_end_index]
	reference_aa_from_codon = CODON_TO_AA.get(codon_sequence, "X")

	if strict_ref_aa and reference_aa_from_codon != protein_change.reference_aa:
		raise click.ClickException(
			"Reference amino acid mismatch: "
			f"HGVS p. requests {protein_change.reference_aa} at position {protein_change.position}, "
			f"but transcript codon {codon_sequence} translates to {reference_aa_from_codon}."
		)

	codon_c_start = codon_cds_start_index + 1

	if is_deletion:
		candidate_variants = [
			CandidateVariant(
				variant_type="del",
				hgvs_c=f"c.{codon_c_start}_{codon_c_start + 2}del",
			)
		]
	else:
		candidate_variants = enumerate_snv_candidates(
			coding_sequence,
			codon_sequence,
			codon_c_start,
			codon_cds_start_index,
			reference_protein,
			protein_change,
		)

		if include_indels:
			candidate_variants.extend(
				enumerate_indel_candidates(
					coding_sequence,
					codon_sequence,
					codon_c_start,
					codon_cds_start_index,
					reference_protein,
					protein_change,
					max_indel_size=max_indel_size,
					use_inv_notation=use_inv_notation,
				)
			)

	deduplicated_candidates: dict[str, CandidateVariant] = {}
	for candidate_variant in candidate_variants:
		deduplicated_candidates[candidate_variant.hgvs_c] = candidate_variant

	sorted_candidates = [deduplicated_candidates[key] for key in sorted(deduplicated_candidates)]

	rows: list[dict[str, str]] = []
	for candidate_variant in sorted_candidates:
		try:
			hgvs_g = map_hgvs_c_to_hgvs_g(
				parser,
				mapper,
				transcript_accession,
				candidate_variant.hgvs_c,
			)
		except Exception as exception:
			click.echo(
				f"Warning: Failed to map {transcript_accession}:{candidate_variant.hgvs_c} to HGVS g. ({exception})",
				err=True,
			)
			hgvs_g = ""

		rows.append(
			{
				"variant_type": candidate_variant.variant_type,
				"hgvs_c": f"{transcript_accession}:{candidate_variant.hgvs_c}",
				"hgvs_g": hgvs_g,
			}
		)

	return rows


def join_variant_rows(
	variant_rows: list[dict[str, str]],
	join_delimiter: str,
) -> dict[str, str]:
	if not variant_rows:
		return {
			"variant_type": "",
			"hgvs_c": "",
			"hgvs_g": "",
		}

	variant_types = [row["variant_type"] for row in variant_rows]
	hgvs_cs = [row["hgvs_c"] for row in variant_rows]
	hgvs_gs = [row["hgvs_g"] for row in variant_rows]

	if len(hgvs_cs) != len(hgvs_gs):
		raise click.ClickException(
			"Internal error: mismatch between hgvs_c and hgvs_g candidate counts for a row."
		)

	return {
		"variant_type": join_delimiter.join(variant_types),
		"hgvs_c": join_delimiter.join(hgvs_cs),
		"hgvs_g": join_delimiter.join(hgvs_gs),
	}


@click.command()
@click.option("--transcript", "transcript_accession", required=False, type=str, help="Coding transcript accession.")
@click.option(
	"--uniprot-id",
	"uniprot_id",
	required=False,
	type=str,
	help="UniProt accession to resolve to MANE Select transcript metadata.",
)
@click.option("--hgvs-p", "hgvs_protein", required=False, type=str, help="Protein HGVS string, e.g. p.Gly12Asp.")
@click.option(
	"--input",
	"input_path",
	type=click.Path(exists=True, dir_okay=False, path_type=Path),
	default=None,
	help="Optional batch input file path. If provided, reads transcript/UniProt/HGVS p. from columns.",
)
@click.option(
	"--input-format",
	default="tsv",
	show_default=True,
	type=click.Choice(["tsv", "csv"], case_sensitive=False),
	help="Batch input file format.",
)
@click.option(
	"--transcript-column",
	default="transcript",
	show_default=True,
	type=str,
	help="Batch mode input column for transcript accession.",
)
@click.option(
	"--uniprot-column",
	default=None,
	type=str,
	help="Batch mode input column for UniProt IDs (used when transcript is absent/blank).",
)
@click.option(
	"--hgvs-p-column",
	default="hgvs_p",
	show_default=True,
	type=str,
	help="Batch mode input column for protein HGVS string.",
)
@click.option(
	"--pass-through-columns",
	"pass_through_selected_columns",
	multiple=True,
	type=str,
	help=(
		"Batch mode: specific additional input column(s) to pass through. "
		"Accepts comma-separated values and may be repeated."
	),
)
@click.option(
	"--pass-through-all-columns",
	is_flag=True,
	default=False,
	help="Batch mode: pass through all additional non-core input columns.",
)
@click.option(
	"--pass-through-prefix",
	default="",
	show_default=True,
	type=str,
	help="Batch mode: prefix for passed-through additional column names.",
)
@click.option(
	"--limit",
	default=None,
	type=click.IntRange(min=1),
	help="Batch mode only: process at most the first N input rows.",
)
@click.option(
	"--uniprot-target",
	default="refseq-protein",
	show_default=True,
	type=click.Choice(["refseq-protein", "ensembl-transcript"], case_sensitive=False),
	help=(
		"UniProt MANE Select identifier to use: RefSeqProteinId (resolved to transcript when possible) "
		"or MANE Ensembl transcript id."
	),
)
@click.option("--assembly", default="GRCh38", show_default=True, type=str, help="Genome assembly for HGVS g. output.")
@click.option(
	"--include-indels",
	is_flag=True,
	default=False,
	help="Include codon-local indels (insertion/deletion/delins) in addition to SNVs.",
)
@click.option(
	"--max-indel-size",
	default=3,
	show_default=True,
	type=click.IntRange(1, 3),
	help="Maximum insertion/deletion size (nt) when --include-indels is enabled.",
)
@click.option(
	"--strict-ref-aa/--no-strict-ref-aa",
	default=True,
	show_default=True,
	help="Require the reference amino acid in HGVS p. to match the transcript codon.",
)
@click.option(
	"--auto-format-hgvs-p",
	is_flag=True,
	default=False,
	help="Auto-convert one-letter amino acid notation (e.g., 'A334D' or 'A334del') to HGVS p. format (e.g., 'p.A334D').",
)
@click.option(
	"--use-inv-notation",
	is_flag=True,
	default=False,
	help="Express eligible reverse-complement delins as HGVS inv variants instead of delins.",
)
@click.option(
	"--one-row-per-input",
	is_flag=True,
	default=False,
	help=(
		"Output one row per input row by joining multiple variant_type/hgvs_c/hgvs_g values "
		"with --join-delimiter."
	),
)
@click.option(
	"--join-delimiter",
	default="|",
	show_default=True,
	type=str,
	help="Delimiter used to join multiple variant_type/hgvs_c/hgvs_g values with --one-row-per-input.",
)
@click.option("--output", type=click.Path(dir_okay=False, path_type=Path), default=None, help="Optional TSV output path.")
def main(
	transcript_accession: str | None,
	uniprot_id: str | None,
	hgvs_protein: str | None,
	input_path: Path | None,
	input_format: str,
	transcript_column: str,
	uniprot_column: str | None,
	hgvs_p_column: str,
	pass_through_selected_columns: tuple[str, ...],
	pass_through_all_columns: bool,
	pass_through_prefix: str,
	limit: int | None,
	uniprot_target: str,
	assembly: str,
	include_indels: bool,
	max_indel_size: int,
	strict_ref_aa: bool,
	auto_format_hgvs_p: bool,
	use_inv_notation: bool,
	one_row_per_input: bool,
	join_delimiter: str,
	output: Path | None,
) -> None:
	uniprot_target = uniprot_target.lower()
	input_format = input_format.lower()

	if not join_delimiter:
		raise click.ClickException("--join-delimiter cannot be empty.")
	if "\t" in join_delimiter:
		raise click.ClickException("--join-delimiter cannot contain a tab character.")
	if "\t" in pass_through_prefix:
		raise click.ClickException("--pass-through-prefix cannot contain a tab character.")
	if pass_through_all_columns and pass_through_selected_columns:
		raise click.ClickException(
			"--pass-through-all-columns cannot be combined with --pass-through-columns. "
			"Use one strategy."
		)
	if pass_through_prefix and not (pass_through_all_columns or pass_through_selected_columns):
		raise click.ClickException(
			"--pass-through-prefix requires --pass-through-all-columns or --pass-through-columns."
		)

	if input_path is None:
		if not hgvs_protein:
			raise click.ClickException("Single mode requires --hgvs-p.")
		if limit is not None:
			raise click.ClickException("--limit is only supported with --input.")

		if auto_format_hgvs_p:
			hgvs_protein = auto_format_hgvs_p_string(hgvs_protein)

		has_transcript = bool((transcript_accession or "").strip())
		has_uniprot = bool((uniprot_id or "").strip())
		if has_transcript == has_uniprot:
			raise click.ClickException(
				"Single mode requires exactly one of --transcript or --uniprot-id, plus --hgvs-p."
			)

	data_provider = hgvs.dataproviders.uta.connect()
	parser = hgvs.parser.Parser()
	mapper = hgvs.assemblymapper.AssemblyMapper(
		data_provider,
		assembly_name=assembly,
		alt_aln_method="splign",
		normalize=True,
	)
	transcript_cache: dict[str, tuple[str, str]] = {}
	uniprot_cache: dict[str, str] = {}

	if input_path is None:
		resolved_transcript_accession = (transcript_accession or "").strip()
		if not resolved_transcript_accession:
			assert uniprot_id is not None
			resolved_transcript_accession = resolve_uniprot_id_to_transcript_accession(
				uniprot_id=uniprot_id,
				uniprot_target=uniprot_target,
				data_provider=data_provider,
				uniprot_cache=uniprot_cache,
			)

		assert hgvs_protein is not None

		rows = reverse_translate_hgvs_p(
			transcript_accession=resolved_transcript_accession,
			hgvs_protein=hgvs_protein,
			include_indels=include_indels,
			max_indel_size=max_indel_size,
			strict_ref_aa=strict_ref_aa,
			use_inv_notation=use_inv_notation,
			parser=parser,
			mapper=mapper,
			data_provider=data_provider,
			transcript_cache=transcript_cache,
		)

		field_names = ["hgvs_p", "transcript", "variant_type", "hgvs_c", "hgvs_g"]
		if one_row_per_input:
			joined_fields = join_variant_rows(rows, join_delimiter=join_delimiter)
			output_rows = [
				{
					"hgvs_p": hgvs_protein,
					"transcript": resolved_transcript_accession,
					"variant_type": joined_fields["variant_type"],
					"hgvs_c": joined_fields["hgvs_c"],
					"hgvs_g": joined_fields["hgvs_g"],
				}
			]
		else:
			output_rows = [
				{
					"hgvs_p": hgvs_protein,
					"transcript": resolved_transcript_accession,
					"variant_type": row["variant_type"],
					"hgvs_c": row["hgvs_c"],
					"hgvs_g": row["hgvs_g"],
				}
				for row in rows
			]

		output_handle = output.open("w", newline="") if output else click.get_text_stream("stdout")
		try:
			writer = csv.DictWriter(output_handle, fieldnames=field_names, delimiter="\t")
			writer.writeheader()
			writer.writerows(output_rows)
		finally:
			if output:
				output_handle.close()

		click.echo(
			f"Generated {len(output_rows)} candidate DNA variants for {hgvs_protein} on {resolved_transcript_accession}.",
			err=True,
		)
		return

	batch_input_delimiter = "\t" if input_format == "tsv" else ","
	input_handle = input_path.open("r", newline="")
	output_handle = output.open("w", newline="") if output else click.get_text_stream("stdout")

	try:
		reader = csv.DictReader(input_handle, delimiter=batch_input_delimiter)
		if reader.fieldnames is None:
			raise click.ClickException(f"Batch input {input_format.upper()} has no header row.")

		has_transcript_column = transcript_column in reader.fieldnames
		has_uniprot_column = bool(uniprot_column) and uniprot_column in reader.fieldnames

		if not has_transcript_column and not has_uniprot_column:
			if not uniprot_column:
				raise click.ClickException(f"Missing transcript column in input {input_format.upper()}: {transcript_column}")
			raise click.ClickException(
				f"Missing transcript and UniProt columns in input {input_format.upper()}: "
				f"{transcript_column}, {uniprot_column}"
			)

		if hgvs_p_column not in reader.fieldnames:
			raise click.ClickException(f"Missing HGVS p column in input {input_format.upper()}: {hgvs_p_column}")

		core_input_columns: list[str] = []
		for column_name in (transcript_column, uniprot_column, hgvs_p_column):
			if (
				column_name
				and column_name in reader.fieldnames
				and column_name not in core_input_columns
			):
				core_input_columns.append(column_name)

		additional_input_columns = [
			column_name for column_name in reader.fieldnames if column_name not in core_input_columns
		]

		selected_passthrough_input_columns: list[str] = []
		if pass_through_all_columns:
			selected_passthrough_input_columns = list(additional_input_columns)
		else:
			seen_selected_passthrough_columns: set[str] = set()
			for raw_column_group in pass_through_selected_columns:
				for raw_column_name in raw_column_group.split(","):
					column_name = raw_column_name.strip()
					if not column_name:
						raise click.ClickException(
							"--pass-through-columns entries cannot be empty. "
							"Use a comma-separated list like 'col1,col2'."
						)
					if column_name in seen_selected_passthrough_columns:
						continue
					seen_selected_passthrough_columns.add(column_name)
					selected_passthrough_input_columns.append(column_name)

			unknown_passthrough_columns = [
				column_name
				for column_name in selected_passthrough_input_columns
				if column_name not in additional_input_columns
			]
			if unknown_passthrough_columns:
				raise click.ClickException(
					"Requested --pass-through-columns not found among additional input columns: "
					+ ", ".join(unknown_passthrough_columns)
				)

		passthrough_output_name_by_input_column: dict[str, str] = {}
		if selected_passthrough_input_columns:
			seen_output_column_names: set[str] = set(core_input_columns)
			seen_output_column_names.update({"variant_type", "hgvs_c", "hgvs_g"})
			for additional_input_column in selected_passthrough_input_columns:
				passthrough_output_column = f"{pass_through_prefix}{additional_input_column}"
				if passthrough_output_column in seen_output_column_names:
					raise click.ClickException(
						"Pass-through column name conflict after applying prefix: "
						f"{additional_input_column} -> {passthrough_output_column}."
					)
				passthrough_output_name_by_input_column[additional_input_column] = passthrough_output_column
				seen_output_column_names.add(passthrough_output_column)

		resolved_transcript_by_uniprot_id: dict[str, str] = {}
		failed_uniprot_ids: set[str] = set()

		field_names = list(core_input_columns)
		if passthrough_output_name_by_input_column:
			field_names.extend(passthrough_output_name_by_input_column.values())
		for column_name in ("variant_type", "hgvs_c", "hgvs_g"):
			if column_name not in field_names:
				field_names.append(column_name)

		writer = csv.DictWriter(output_handle, fieldnames=field_names, delimiter="\t")
		writer.writeheader()

		total_input_rows = 0
		total_output_rows = 0
		total_rows_with_variants = 0
		total_rows_with_errors = 0
		total_rows_skipped_missing_hgvs_p = 0

		for row_index, row in enumerate(reader, start=1):
			if limit is not None and total_input_rows >= limit:
				break

			total_input_rows += 1
			row_output_template = {column_name: (row.get(column_name) or "") for column_name in core_input_columns}
			if passthrough_output_name_by_input_column:
				for input_column_name, output_column_name in passthrough_output_name_by_input_column.items():
					row_output_template[output_column_name] = row.get(input_column_name) or ""

			row_transcript = (row.get(transcript_column) or "").strip() if has_transcript_column else ""
			row_uniprot_id = (row.get(uniprot_column) or "").strip() if has_uniprot_column and uniprot_column else ""

			if not row_transcript and row_uniprot_id:
				if row_uniprot_id in resolved_transcript_by_uniprot_id:
					row_transcript = resolved_transcript_by_uniprot_id[row_uniprot_id]
				elif row_uniprot_id not in failed_uniprot_ids:
					try:
						row_transcript = resolve_uniprot_id_to_transcript_accession(
							uniprot_id=row_uniprot_id,
							uniprot_target=uniprot_target,
							data_provider=data_provider,
							uniprot_cache=uniprot_cache,
						)
						resolved_transcript_by_uniprot_id[row_uniprot_id] = row_transcript
					except Exception as exception:
						failed_uniprot_ids.add(row_uniprot_id)
						click.echo(
							f"Warning: Failed to resolve UniProt ID {row_uniprot_id} ({exception}).",
							err=True,
						)
				if has_transcript_column:
					row_output_template[transcript_column] = row_transcript

			row_hgvs_p = (row.get(hgvs_p_column) or "").strip()
			if auto_format_hgvs_p:
				row_hgvs_p = auto_format_hgvs_p_string(row_hgvs_p)

			if not row_hgvs_p:
				total_rows_skipped_missing_hgvs_p += 1
				continue

			if not row_transcript:
				click.echo(
					f"Warning: Row {row_index} is missing transcript/UniProt-derived transcript; "
					"writing empty variant fields.",
					err=True,
				)
				total_rows_with_errors += 1
				writer.writerow({**row_output_template, "variant_type": "", "hgvs_c": "", "hgvs_g": ""})
				total_output_rows += 1
				continue

			try:
				row_variants = reverse_translate_hgvs_p(
					transcript_accession=row_transcript,
					hgvs_protein=row_hgvs_p,
					include_indels=include_indels,
					max_indel_size=max_indel_size,
					strict_ref_aa=strict_ref_aa,
					use_inv_notation=use_inv_notation,
					parser=parser,
					mapper=mapper,
					data_provider=data_provider,
					transcript_cache=transcript_cache,
				)
			except Exception as exception:
				click.echo(
					f"Warning: Row {row_index} failed reverse translation ({exception}); writing empty variant fields.",
					err=True,
				)
				total_rows_with_errors += 1
				writer.writerow({**row_output_template, "variant_type": "", "hgvs_c": "", "hgvs_g": ""})
				total_output_rows += 1
				continue

			if row_variants:
				total_rows_with_variants += 1
				if one_row_per_input:
					joined_fields = join_variant_rows(row_variants, join_delimiter=join_delimiter)
					writer.writerow(
						{
							**row_output_template,
							"variant_type": joined_fields["variant_type"],
							"hgvs_c": joined_fields["hgvs_c"],
							"hgvs_g": joined_fields["hgvs_g"],
						}
					)
					total_output_rows += 1
				else:
					for row_variant in row_variants:
						writer.writerow(
							{
								**row_output_template,
								"variant_type": row_variant["variant_type"],
								"hgvs_c": row_variant["hgvs_c"],
								"hgvs_g": row_variant["hgvs_g"],
							}
						)
						total_output_rows += 1
			else:
				writer.writerow({**row_output_template, "variant_type": "", "hgvs_c": "", "hgvs_g": ""})
				total_output_rows += 1

		if total_rows_skipped_missing_hgvs_p > 0:
			click.echo(
				"Warning: Skipped "
				f"{total_rows_skipped_missing_hgvs_p} row(s) with missing HGVS p. values.",
				err=True,
			)

		click.echo(
			"Batch reverse-translation summary: "
			f"{total_input_rows} input rows, "
			f"{total_output_rows} output rows, "
			f"{total_rows_with_variants} rows with >=1 variant, "
			f"{total_rows_with_errors} rows with warnings/errors, "
			f"{total_rows_skipped_missing_hgvs_p} rows skipped for missing HGVS p..",
			err=True,
		)
	finally:
		input_handle.close()
		if output:
			output_handle.close()


if __name__ == "__main__":
	main()
