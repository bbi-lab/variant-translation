# variant-translation: Tools for translating representations of genetic variants

## Overview

`variant-translation` provides Python and CLI workflows for translating and comparing
HGVS variant representations in large datasets.

At a high level, this repository includes:

- **Reverse translation**: convert protein HGVS changes (`hgvs_p`) into candidate
   coding and genomic variants (`hgvs_c`, `hgvs_g`) using transcript-aware mapping.
- **File comparison**: compare two reverse-translated result files and report rows
   unique to each file plus rows with differing variant outputs.
- **Batch utilities**: stream-friendly TSV/CSV processing features (skip/limit,
   pass-through columns, one-row-per-input output, error-row capture) for
   production-scale pipelines.

The packaged CLI entry points are:

- `reverse-translate-variants`
- `compare-reverse-translated-variants`

The repository also contains small standalone helper scripts under `scripts/` for
one-off data wrangling tasks.

## Quick start

After installing dependencies, run one of these examples:

- Reverse translate a single protein variant:

   ```bash
   python -m src.scripts.reverse_translate_variants --transcript NM_000546.6 --hgvs-p p.Arg175His
   ```

- Reverse translate a batch TSV file:

   ```bash
   python -m src.scripts.reverse_translate_variants \
       --input input.tsv \
       --input-format tsv \
       --hgvs-p-column hgvs_p \
       --transcript-column transcript \
       --output reverse_translated.tsv \
       --errors reverse_translation_errors.tsv
   ```

- Compare two reverse-translated files:

   ```bash
   python -m src.scripts.compare_reverse_translated_variants \
       --a-input a.tsv \
       --b-input b.tsv \
       --a-only-output a_only.tsv \
       --b-only-output b_only.tsv \
       --different-output differences.tsv
   ```

## Installation

### Prerequisites

- Python 3.8 or later
- Virtual environment (recommended)

### Option A: Install from PyPI (library + CLI)

1. Create and activate a virtual environment:

   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On macOS/Linux
   # or
   .venv\Scripts\activate  # On Windows
   ```

2. Install from PyPI:

   ```bash
   pip install variant-translation
   ```

3. Run installed CLI entry points:

   ```bash
   reverse-translate-variants --help
   compare-reverse-translated-variants --help
   ```

### Option B: Install from GitHub source

1. Clone this repository:

   ```bash
   git clone git@github.com:bbi-lab/variant-translation.git
   cd variant-translation
   ```

2. Create and activate a virtual environment:

   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On macOS/Linux
   # or
   .venv\Scripts\activate  # On Windows
   ```

3. Install dependencies:

   ```bash
   pip install -r requirements.txt
   ```

4. (Optional) Install as an editable package to use CLI entry points while developing:

   ```bash
   pip install -e .
   ```

## Reverse-Translate Protein Variants

### `reverse_translate_variants.py`

Enumerate DNA variants that produce a requested protein HGVS change:

If installed from PyPI, use the `reverse-translate-variants` command in place of
`python -m src.scripts.reverse_translate_variants`.

#### Programmatic usage (without calling `main`)

You can call the core reverse-translation function directly from Python and reuse initialized HGVS/UTA objects
across multiple variants:

```python
import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.parser

from src.scripts.reverse_translate_variants import (
   auto_format_hgvs_p_string,
   join_variant_rows,
   reverse_translate_batch_rows,
   reverse_translate_hgvs_p,
)

data_provider = hgvs.dataproviders.uta.connect()
parser = hgvs.parser.Parser()
mapper = hgvs.assemblymapper.AssemblyMapper(
   data_provider,
   assembly_name="GRCh38",
   alt_aln_method="splign",
   normalize=True,
)

transcript_cache: dict[str, tuple[str, str]] = {}

hgvs_p = auto_format_hgvs_p_string("R175H")
rows = reverse_translate_hgvs_p(
   transcript_accession="NM_000546.6",
   hgvs_protein=hgvs_p,
   include_indels=True,
   max_indel_size=3,
   strict_ref_aa=True,
   use_inv_notation=False,
   allow_length_changing_stop_candidates=True,
   parser=parser,
   mapper=mapper,
   data_provider=data_provider,
   transcript_cache=transcript_cache,
)

print(rows[0])
# {'variant_type': 'snv', 'hgvs_c': 'NM_000546.6:c....', 'hgvs_g': 'NC_...:g....'}

# Optional: produce one-row-per-input style joined fields
joined = join_variant_rows(rows, join_delimiter="|")
print(joined)
```

#### Programmatic batch loop (without calling CLI `main`)

For batch processing, use `reverse_translate_batch_rows(...)` to process row iterables directly. Initialize
UTA/HGVS once and reuse `parser`, `mapper`, and `transcript_cache` for every row:

```python
import csv

import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.parser

from src.scripts.reverse_translate_variants import (
   reverse_translate_batch_rows,
)

data_provider = hgvs.dataproviders.uta.connect()
parser = hgvs.parser.Parser()
mapper = hgvs.assemblymapper.AssemblyMapper(
   data_provider,
   assembly_name="GRCh38",
   alt_aln_method="splign",
   normalize=True,
)

transcript_cache: dict[str, tuple[str, str]] = {}

with open("input.csv", newline="") as in_handle:
    reader = csv.DictReader(in_handle)
    output_rows = reverse_translate_batch_rows(
        rows=reader,
        transcript_column="Ref_seq_transcript_ID",
        hgvs_p_column="hgvs_p",
        include_indels=True,
        max_indel_size=3,
        strict_ref_aa=True,
        use_inv_notation=True,
        allow_length_changing_stop_candidates=False,
        parser=parser,
        mapper=mapper,
        data_provider=data_provider,
        transcript_cache=transcript_cache,
        auto_format_hgvs_p=True,
        one_row_per_input=True,
        join_delimiter="^",
        skip_missing_hgvs_p=False,
    )

with open("output.tsv", "w", newline="") as out_handle:
    output_columns = [
        "ID",
        "Dataset",
        "Gene",
        "Ref_seq_transcript_ID",
        "hgvs_p",
        "variant_type",
        "hgvs_c",
        "hgvs_g",
    ]
    writer = csv.DictWriter(out_handle, fieldnames=output_columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(output_rows)
```

`reverse_translate_batch_rows(...)` options mirror key CLI behavior:
- `one_row_per_input` + `join_delimiter`: same output shape as `--one-row-per-input` / `--join-delimiter`
- `auto_format_hgvs_p`: same behavior as `--auto-format-hgvs-p`
- `skip_missing_hgvs_p`: default `True` (CLI-like skip behavior), set `False` to emit blank variant fields
- `raise_on_error`: default `False`; set `True` for fail-fast pipelines that should stop on the first
   reverse-translation error
- `error_rows`: optional list to collect rows that encountered errors; each entry includes the same
   output columns plus an `error` message

Useful functions:
- `reverse_translate_hgvs_p(...)`: core reverse-translation function returning candidate rows with
   `variant_type`, `hgvs_c`, and `hgvs_g`
- `reverse_translate_batch_rows(...)`: batch helper for iterables of row dictionaries
- `auto_format_hgvs_p_string(...)`: converts shorthand protein strings (e.g., `R175H`, `A334del`) to HGVS p. form
- `join_variant_rows(...)`: joins multi-candidate output into single delimited fields (same behavior as
   `--one-row-per-input`)

> **Requires UTA.** This script connects to a UTA (Universal Transcript Archive) PostgreSQL database to look up
> transcript sequences and map c. variants to g. variants. See [UTA Setup](#uta-universal-transcript-archive) in
> the Configuration section.

```bash
python -m src.scripts.reverse_translate_variants --transcript NM_000546.6 --hgvs-p p.Arg175His
```

You can also provide a UniProt ID instead of a transcript accession. The script resolves the UniProt entry to a
MANE Select identifier before reverse translation:

```bash
python -m src.scripts.reverse_translate_variants --uniprot-id P04637 --hgvs-p p.Arg175His
```

By default, UniProt resolution uses `RefSeqProteinId` from
`uniProtKBCrossReferences[0].properties` and resolves it to a transcript accession. To use the MANE Select Ensembl
transcript directly (`uniProtKBCrossReferences[0].id`), set:

```bash
python -m src.scripts.reverse_translate_variants --uniprot-id P04637 --hgvs-p p.Arg175His --uniprot-target ensembl-transcript
```

**How transcript selection works**

- **Single mode (`--input` not provided):** you must supply `--hgvs-p` and **exactly one** of
   `--transcript` or `--uniprot-id`.
- **Batch mode (`--input` provided):** the file must include `--hgvs-p-column` and at least one transcript source:
   `--transcript-column` (default: `transcript`) and/or `--uniprot-column`.
- If both transcript and UniProt columns exist, each row uses the transcript value when present; UniProt resolution is
   used as a fallback when transcript is blank.
- Rows with neither a transcript nor a resolvable UniProt-derived transcript are written with empty variant fields and
   a warning (and can be captured in `--errors`).

**Common transcript/UniProt mode errors**

- `Single mode requires --hgvs-p.`
   - Add `--hgvs-p p.<change>` when running without `--input`.
- `Single mode requires exactly one of --transcript or --uniprot-id, plus --hgvs-p.`
   - In single mode, pass only one transcript source: either `--transcript` or `--uniprot-id`.
- `Missing transcript column in input ...` or `Missing transcript and UniProt columns in input ...`
   - In batch mode, ensure your header contains the column named by `--transcript-column` and/or
     `--uniprot-column`.
- `Missing HGVS p column in input ...`
   - Ensure your input header contains the column named by `--hgvs-p-column`.
- `--limit is only supported with --input.` / `--skip is only supported with --input.` / `--errors is only supported with --input.`
   - These options are batch-mode only. Add `--input ...` or remove those flags in single-variant mode.

By default, this returns all codon-local SNVs that yield the same amino-acid change, along with HGVS c. and HGVS g.
expressions.

To also consider small indels (up to 3 nt by default), including full-codon replacements
(`c.N_N+2delinsNNN`) for substitutions that are not SNV-accessible:

```bash
python -m src.scripts.reverse_translate_variants --transcript NM_000546.6 --hgvs-p p.Arg175His --include-indels
```

To express eligible reverse-complement replacements as HGVS inversions (`inv`) instead of `delins`:

```bash
python -m src.scripts.reverse_translate_variants --transcript NM_000546.6 --hgvs-p p.Arg175His --include-indels --use-inv-notation
```

To restrict premature-stop reverse translation to substitutions and same-length `delins` variants, excluding
length-changing insertion/deletion candidates:

```bash
python -m src.scripts.reverse_translate_variants --transcript NM_000546.6 --hgvs-p p.Arg175Ter --include-indels --substitutions-and-delins-only
```

**Representation notes**

- For downstream use, it is still a good idea to pass generated variants through
   [VariantValidator](https://variantvalidator.org/) or a similar HGVS-aware normalization/validation tool to confirm
   canonical form. In current testing we have not observed generation of less-recommended representations, but
   external validation remains recommended for production pipelines.
- During c.→g. projection, if the affected codon spans an intron, the genomic expression may include an intronic
   deletion segment. This can still be a valid reverse translation of the same coding effect. In these cases there are
   often many genomic alternatives that reduce to the same c. variant, and it may be preferable to represent a
   multi-variant/compound genomic change that alters only the coding bases on each side of the exon boundary.
   - Toy pattern: a boundary-spanning genomic edit can appear as a single contiguous event
     (`NC_...:g.<left>_<right>delins...`), while an alternative representation with the same coding consequence may be
     a compound form (`NC_...:g.[<left_exonic_change>;<right_exonic_change>]`) that avoids explicitly deleting intronic
     sequence.

**Supported protein HGVS inputs**

This tool currently accepts **single-amino-acid changes at one position**:

- missense/substitution, e.g. `p.Arg175His`
- synonymous, e.g. `p.Arg175=`
- premature stop / stop-gain, e.g. `p.Arg175Ter` or `p.Arg175*`
- single-amino-acid deletion, e.g. `p.Arg175del`

It also accepts the same forms with an optional protein accession prefix, e.g.
`NP_000537.3:p.Arg175His`, and `--auto-format-hgvs-p` can help with shorthand
inputs such as `R175H`, `R175=`, `A334del`, or `A334-`.

When `--include-indels` is enabled, the script may return DNA-level
insertions, deletions, `delins`, or `inv` candidates **that produce one of the
supported protein outcomes above**. That option does **not** mean the input
protein HGVS can itself be an arbitrary protein indel notation.

**Not currently supported**

The reverse-translation parser does not currently handle broader protein HGVS
classes such as:

- single-amino-acid insertions
- multi-residue deletions, insertions, duplications, or `delins`
- frameshifts (for example `fs` forms)
- extensions / stop-loss forms (for example `ext*?`); note that plain stop-loss
   requests with `*` as the reference amino acid are explicitly rejected
- complex protein expressions spanning more than one amino-acid position

For those cases, the script will either reject the input as unparseable or, for
stop-loss requests, report that the variant type is not supported.

**Key options:**
- `--transcript`: Single mode transcript accession (required in single mode unless using `--uniprot-id`)
- `--transcript-column`: Batch mode input column containing transcript accessions (default: `transcript`)
- `--assembly`: Genome assembly for HGVS g. projection (default: `GRCh38`)
- `--uniprot-id`: Resolve UniProt accession to a MANE Select identifier for single-variant mode
- `--uniprot-column`: Batch mode input column for UniProt IDs used when transcript is blank/missing
- `--uniprot-target`: UniProt MANE identifier source: `refseq-protein` (default) or `ensembl-transcript`
- `--include-indels`: Include codon-local insertion/deletion/delins candidates
- `--use-inv-notation`: Express eligible reverse-complement delins as HGVS `inv` variants
- `--substitutions-and-delins-only`: For premature stops, suppress length-changing insertion/deletion candidates
- `--skip`: Batch mode only; skip the first N input rows after the header
- `--max-indel-size`: Maximum insertion/deletion size in nucleotides (1-3)
- `--limit`: Batch mode only; process at most the first N input rows
- `--input-format`: Batch input format: `tsv` (default) or `csv`
- `--csv-field-size-limit`: Batch mode parser field-length limit (default: current Python `csv` module limit)
- `--pass-through-all-columns`: Batch mode; pass through all additional non-core input columns
- `--pass-through-columns`: Batch mode; pass through specific additional column(s) (comma-separated and repeatable)
- `--pass-through-prefix`: Batch mode; prefix for passed-through additional column names
- `--one-row-per-input`: Emit exactly one output row per input row by joining multiple candidates
- `--join-delimiter`: Delimiter used with `--one-row-per-input` (default: `|`)
- `--output`: Write TSV output to a file (otherwise prints to stdout)
- `--errors`: Batch mode only; write warning/error rows to a TSV file with an added `error` column

**Batch mode (default TSV):**

```bash
python -m src.scripts.reverse_translate_variants --input input.tsv --output output.tsv
```

By default, batch mode reads transcript accessions from `transcript` and protein HGVS from `hgvs_p` (override with
`--transcript-column` and `--hgvs-p-column`).

If your file does not have transcript accessions, provide `--uniprot-column` so transcript IDs can be resolved from
UniProt. At least one of `--transcript-column` or `--uniprot-column` must be present in the input header.

The output includes core columns (`transcript`, optional `uniprot`, and `hgvs_p`) plus `variant_type`, `hgvs_c`, and
`hgvs_g`.

Rows with missing/blank `hgvs_p` values are skipped. If one or more rows are skipped, the script emits a warning.

To pass through additional input columns, enable:

```bash
python -m src.scripts.reverse_translate_variants --input input.tsv --pass-through-all-columns --output output.tsv
```

To pass through only selected additional columns:

```bash
python -m src.scripts.reverse_translate_variants \
   --input input.tsv \
   --pass-through-columns sample_id,project_id \
   --pass-through-columns plate_id \
   --output output.tsv
```

To prefix passed-through additional column names (works with either strategy):

```bash
python -m src.scripts.reverse_translate_variants --input input.tsv --pass-through-all-columns --pass-through-prefix input_ --output output.tsv
```

To drive batch mode from UniProt IDs instead of transcript accessions:

```bash
python -m src.scripts.reverse_translate_variants \
   --input input.tsv \
   --uniprot-column uniprot_id \
   --output output.tsv
```

To process only the first N rows in batch mode (after the header):

```bash
python -m src.scripts.reverse_translate_variants --input input.tsv --limit 100 --output output.tsv
```

To skip the first N rows in batch mode (after the header), use `--skip N`. For example, `--skip 100` skips input
lines 1 through 101 (header + first 100 records):

```bash
python -m src.scripts.reverse_translate_variants --input input.tsv --skip 100 --output output.tsv
```

To write warning/error rows (for example, failed reverse translation rows) to a separate file with an
`error` message column:

```bash
python -m src.scripts.reverse_translate_variants --input input.tsv --output output.tsv --errors errors.tsv
```

If your batch input has unusually long fields, increase the parser limit:

```bash
python -m src.scripts.reverse_translate_variants --input input.tsv --csv-field-size-limit 1000000 --output output.tsv
```

In UniProt batch mode, the script processes input rows as a stream and resolves UniProt IDs on demand, caching
resolved transcripts so each unique ID is resolved at most once.

To output exactly one row per input row (joining multiple candidate values in each of `variant_type`, `hgvs_c`, and
`hgvs_g`):

```bash
python -m src.scripts.reverse_translate_variants \
   --input input.tsv \
   --uniprot-column uniprot_id \
   --hgvs-p-column hgvs_p \
   --include-indels \
   --one-row-per-input \
   --join-delimiter "|" \
   --output output_joined.tsv
```

To use CSV input instead of TSV:

```bash
python -m src.scripts.reverse_translate_variants --input input.csv --input-format csv --output output.tsv
```

In this mode, the joined `hgvs_c` and `hgvs_g` lists are kept aligned so each position corresponds to the same
candidate variant.

**Auto-formatting HGVS p. strings:**

If your HGVS p. column contains one-letter amino acid notation (e.g., `A334D`, `M1V`) or single-amino-acid deletions
(e.g., `A334del`, `A334-`), use the `--auto-format-hgvs-p` flag to automatically convert these to HGVS p. format
(e.g., `p.A334D`, `p.M1V`, `p.A334del`). This applies to both single-mode and batch-mode:

```bash
# Single mode: auto-format A334D → p.A334D before processing
python -m src.scripts.reverse_translate_variants --transcript NM_000546.6 --hgvs-p A334D --auto-format-hgvs-p

# Batch mode: auto-format one-letter notations in the hgvs_p column
python -m src.scripts.reverse_translate_variants --input input.tsv --auto-format-hgvs-p --output output.tsv
```

Supported patterns: one-letter (e.g., `A334D`, `M1V`) and three-letter (e.g., `Ala334Asp`, `Met1Val`) amino acid
codes; single stop codon changes (e.g., `*123C`); synonymous changes (e.g., `R175=`); and deletions (e.g., `A334del`,
`A334-`, `Met1del`).

### Compare Reverse-Translated Variant Files

Compare two single-line reverse-translation files (A and B), matching rows by exact transcript + `hgvs_p` and
reporting A-only keys, B-only keys, and rows where reverse translations differ.

If installed from PyPI, use the `compare-reverse-translated-variants` command in place of
`python -m src.scripts.compare_reverse_translated_variants`.

```bash
python -m src.scripts.compare_reverse_translated_variants \
   --a-input a.tsv \
   --b-input b.tsv \
   --a-only-output a_only.tsv \
   --b-only-output b_only.tsv \
   --different-output differences.tsv \
   --a-missing-output a_missing.tsv \
   --b-missing-output b_missing.tsv
```

Key options:
- Per-file formats: `--a-input-format`, `--b-input-format` (`tsv`/`csv`)
- Per-file single-line separators: `--a-separator`, `--b-separator`
- Per-file column options: transcript, `hgvs_p`, `hgvs_c`, `hgvs_g`
- Additional key columns for matching: `--match-columns` (comma-separated and repeatable)
- Comparison scope: `--compare hgvs_c|hgvs_g|both`
- Reference-insensitive comparison: `--ignore-reference-differences`
- Indel-normalized comparison: `--normalize-indels` to normalize HGVS variants before comparison so equivalent indel/codon-replacement representations can match
- Differences output pass-through:
   - `--pass-through-all-columns` for all non-core columns
   - `--pass-through-columns` for selected non-core columns from both A and B (comma-separated and repeatable)
   - `--a-pass-through-columns` and `--b-pass-through-columns` for input-specific selected non-core columns
   - `--a-pass-through-prefix` and `--b-pass-through-prefix` for output naming
- Missing-key outputs: `--a-missing-output` and `--b-missing-output` for rows missing transcript or `hgvs_p`

If one file represents an event as a short indel and the other represents the same biological change as a
full-codon `delins`, use `--normalize-indels` to normalize each HGVS expression through the `hgvs` library before
comparison. This is slower than the default exact-string comparison and requires access to UTA (and benefits from a
local SeqRepo), so it is only enabled when explicitly requested.

Example using separate pass-through selections for each input:

```bash
python -m src.scripts.compare_reverse_translated_variants \
   --a-input a.tsv \
   --b-input b.tsv \
   --a-pass-through-columns sample_id,plate_id \
   --b-pass-through-columns record_id,batch \
   --a-only-output a_only.tsv \
   --b-only-output b_only.tsv \
   --different-output differences.tsv
```

Example allowing equivalent indel representations to match:

```bash
python -m src.scripts.compare_reverse_translated_variants \
   --a-input a.tsv \
   --b-input b.tsv \
   --compare hgvs_c \
   --ignore-reference-differences \
   --normalize-indels \
   --a-only-output a_only.tsv \
   --b-only-output b_only.tsv \
   --different-output differences.tsv
```

## Configuration

### UTA (Universal Transcript Archive)

The `reverse_translate_variants` script depends on a UTA PostgreSQL database, which provides transcript sequences
and genomic coordinate mappings for the [biocommons `hgvs`](https://github.com/biocommons/hgvs) library. By default
the `hgvs` library connects to the public biocommons cloud instance, which can be slow or rate-limited. Running a
local Docker container is strongly recommended for any non-trivial use.

#### Step 1 — Pull and start the UTA Docker container

The `biocommons/uta` image provides a PostgreSQL instance but ships with an **empty database**. Data must be loaded
separately from a dump file (see Step 3).

```bash
docker pull biocommons/uta:uta_20241220

docker run \
  --name uta \
  --platform linux/amd64 \
  -e POSTGRES_PASSWORD=uta_admin \
  -p 5433:5432 \
  -d biocommons/uta:uta_20241220
```

Port `5433` is used here to avoid conflicts with any other local PostgreSQL instance on the default port `5432`.
Change it to any free port you prefer.

Wait for PostgreSQL to finish initialising before proceeding. Check the logs until you see
`database system is ready to accept connections`:

```bash
docker logs -f uta
```

#### Step 2 — Download the UTA database dump

The dump file is available from the biocommons download server (~1–2 GB). Download it with `curl` or `wget`:

```bash
curl -LO "http://dl.biocommons.org/uta/uta_20241220.pgd.gz"
```

If the download fails or is very slow, you can also find the file on the
[biocommons/uta GitHub releases page](https://github.com/biocommons/uta/releases).

#### Step 3 — Restore the dump into the container

The dump is a gzip-compressed plain-SQL file. Pipe it directly into the running container:

```bash
gzip -cdq uta_20241220.pgd.gz | docker exec -i uta psql -U uta_admin uta
```

This will take several minutes to complete. There is no progress indicator; the shell prompt will return when the
restore is done.

#### Step 4 — Set the UTA_DB_URL environment variable

Tell the `hgvs` library (and therefore this script) where to find the local instance:

```bash
export UTA_DB_URL="postgresql://anonymous@localhost:5433/uta/uta_20241220"
```

You can add this line to your shell profile or `.env` file so it persists across sessions.

If `UTA_DB_URL` is not set the library falls back to the public cloud instance at
`postgresql://anonymous:anonymous@uta.biocommons.org/uta/uta_20241220`.

#### Step 5 — Optional: use a local SeqRepo to avoid repeated network fetches

When mapping c. variants to g. coordinates, the `hgvs` library calls
`SeqFetcher.fetch_seq` for each candidate variant. Without a local sequence
store this makes a fresh NCBI/Ensembl HTTP request for every call — there is
**no internal caching** in `SeqFetcher` or in the UTA data provider's
`get_seq` method.

> **Note:** The script's own `transcript_cache` dict prevents repeated
> `get_seq` calls for the same transcript accession, but it does not cover the
> genomic-sequence lookups made internally by `AssemblyMapper` during c→g
> mapping. Those are what benefit most from a local SeqRepo.

Install and populate a local [SeqRepo](https://github.com/biocommons/biocommons.seqrepo)
snapshot, then set `HGVS_SEQREPO_DIR` so the `hgvs` library uses it instead
of the network:

```bash
# Install the seqrepo command-line tool
pip install seqrepo

# Download the latest snapshot (several GB; run once)
sudo seqrepo --root-directory /usr/local/share/seqrepo pull

# The pull creates a dated directory and a 'latest' symlink, e.g.:
#   /usr/local/share/seqrepo/2021-01-29/   ← snapshot
#   /usr/local/share/seqrepo/latest        ← symlink → 2021-01-29
export HGVS_SEQREPO_DIR=/usr/local/share/seqrepo/latest
```

You can add this variable to your shell profile or `.env` file so it persists:

```dotenv
export HGVS_SEQREPO_DIR=/usr/local/share/seqrepo/latest
```

With `HGVS_SEQREPO_DIR` set, the `hgvs` `SeqFetcher` switches to a local
`biocommons.seqrepo.SeqRepo` instance and serves all sequences from disk
without any network calls, which significantly speeds up large batch runs.

#### Stop, restart, and remove the container

```bash
docker stop uta
docker start uta   # reuses the existing container; no need to restore again
```

To remove the container entirely (data will be lost; restore step must be repeated):

```bash
docker rm -f uta
```

#### Verify the connection

```bash
psql "postgresql://anonymous@localhost:5433/uta" -c "SELECT count(*) FROM uta_20241220.tx_exon_aln_v LIMIT 1;"
```

A row count result confirms the database is accessible and the schema is populated.

## Dependencies

Key dependencies include:

- **click**: Command-line interface creation
- **hgvs**: HGVS variant parsing and validation

See [requirements.txt](requirements.txt) for the complete list.

## Publish to PyPI

From a clean, version-tagged working tree:

Before building, bump `version` in `pyproject.toml`.

### Release checklist

- [ ] Bump `version` in `pyproject.toml`.
- [ ] Ensure tests pass (`pytest`).
- [ ] Rebuild distributions from a clean state (`rm -rf dist build *.egg-info`).
- [ ] Build artifacts (`python -m build`).
- [ ] Validate metadata (`python -m twine check dist/*`).
- [ ] (Recommended) Upload and test install via TestPyPI.
- [ ] Upload to PyPI (`python -m twine upload dist/*`).

### One-command release targets

This repository includes a `Makefile` with release targets that implement the checklist steps:

```bash
make release-prep        # test + clean + build + twine check
make release-testpypi    # release-prep + upload to TestPyPI
make release             # release-prep + upload to PyPI
```

1. Install publishing dependencies:

   ```bash
   pip install .[publish]
   ```

2. Build source and wheel distributions:

   ```bash
   python -m build
   ```

3. Validate distribution metadata:

   ```bash
   python -m twine check dist/*
   ```

4. Upload to PyPI:

   ```bash
   python -m twine upload dist/*
   ```

### TestPyPI dry-run (recommended before production upload)

Use TestPyPI to validate packaging and installation without publishing to the main index:

1. Upload to TestPyPI:

   ```bash
   python -m twine upload --repository testpypi dist/*
   ```

2. Install from TestPyPI in a clean virtual environment:

   ```bash
   pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple variant-translation
   ```
