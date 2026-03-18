# variant-translation: Tools for translating representations of genetic variants

## Installation

### Prerequisites

- Python 3.8 or later
- Virtual environment (recommended)

### Setup

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

## Reverse-Translate Protein Variants

### `reverse_translate_variants.py`

Enumerate DNA variants that produce a requested protein HGVS change:

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

By default, this returns all codon-local SNVs that yield the same amino-acid change, along with HGVS c. and HGVS g.
expressions.

To also consider small indels (up to 3 nt by default), including full-codon replacements
(`c.N_N+2delinsNNN`) for substitutions that are not SNV-accessible:

```bash
python -m src.scripts.reverse_translate_variants --transcript NM_000546.6 --hgvs-p p.Arg175His --include-indels
```

**Key options:**
- `--assembly`: Genome assembly for HGVS g. projection (default: `GRCh38`)
- `--uniprot-id`: Resolve UniProt accession to a MANE Select identifier for single-variant mode
- `--uniprot-target`: UniProt MANE identifier source: `refseq-protein` (default) or `ensembl-transcript`
- `--include-indels`: Include codon-local insertion/deletion/delins candidates
- `--max-indel-size`: Maximum insertion/deletion size in nucleotides (1-3)
- `--limit`: Batch mode only; process at most the first N input rows
- `--input-format`: Batch input format: `tsv` (default) or `csv`
- `--pass-through-all-columns`: Batch mode; pass through all additional non-core input columns
- `--pass-through-columns`: Batch mode; pass through specific additional column(s) (comma-separated and repeatable)
- `--pass-through-prefix`: Batch mode; prefix for passed-through additional column names
- `--one-row-per-input`: Emit exactly one output row per input row by joining multiple candidates
- `--join-delimiter`: Delimiter used with `--one-row-per-input` (default: `|`)
- `--output`: Write TSV output to a file (otherwise prints to stdout)

**Batch mode (default TSV):**

```bash
python -m src.scripts.reverse_translate_variants --input input.tsv --output output.tsv
```

By default, batch mode expects columns named `transcript` and `hgvs_p` (override with
`--transcript-column` and `--hgvs-p-column`). The output includes core columns (`transcript`, optional `uniprot`, and
`hgvs_p`) plus `variant_type`, `hgvs_c`, and `hgvs_g`.

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
- Differences output pass-through:
   - `--pass-through-all-columns` for all non-core columns
   - `--pass-through-columns` for selected non-core columns (comma-separated and repeatable)
   - `--a-pass-through-prefix` and `--b-pass-through-prefix` for output naming
- Missing-key outputs: `--a-missing-output` and `--b-missing-output` for rows missing transcript or `hgvs_p`

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
