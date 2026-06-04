# Pooled Ordering from Raw Sequences Design

**Date:** 2026-06-03
**Branch:** `parker/update-infra` (or new branch off master)
**Scope:** New `pyguide-collate-seq` command + `pooled-seq` order format

---

## 1. Goals

Allow users to supply raw 20nt sgRNA spacer sequences directly (instead of gene names or Horlbeck library guide IDs) and produce a pooled Agilent oligo order file with restriction sites and primers — using the same downstream pipeline as the existing `pooled` order format.

---

## 2. User Workflow

```bash
# Step 1: collate sequence file(s) with primer assignment
pyguide-collate-seq --sequence_files my_guides_lib1.txt my_guides_lib2.txt

# Step 2: order (mirrors existing pooled flow, new format flag)
pyguide-order --wishlist_file collated_seq_wishlist_<date>.txt \
              --name Name --ai i --order_format pooled-seq
```

---

## 3. Components

### 3.1 New file: `pyguide/collate_seq.py`

Contains all logic for the `pyguide-collate-seq` command.

**Functions:**

`validate_sequence_file(file: str) -> pd.DataFrame`
- Reads a tab-delimited file (no header): column 1 = name, column 2 = sequence
- Collects ALL validation errors before raising, so user sees every problem at once
- Validation rules:
  - Sequence length must be exactly 20nt
  - Sequence must contain only A/C/G/T characters (case-insensitive; normalizes to uppercase)
  - No duplicate names within the file
  - File must not be empty
- On error: prints all invalid rows with row number, name, sequence, and specific problem; exits with code 1
- On success: returns a DataFrame with columns `name`, `seq`

`generate_pooled_seq_list(sequence_files: List[str], primer_list: List[Tuple[str, str]]) -> None`
- For each (sequence_file, primer_tuple) pair:
  - Calls `validate_sequence_file` to read and validate
  - Writes rows as `name\tseq\tleft_primer\tright_primer\tlib_num` to the collated output file
- Output filename: `collated_seq_wishlist_<YYMMDD>.txt` in the same directory as the first input file
- Uses `get_unique_filename` from `guide.py` to avoid overwriting existing files
- Reuses `get_primer_list` from `pool.py` for primer assignment

`main(raw_args=None)`
- CLI entry point for `pyguide-collate-seq`
- Arguments:
  - `--sequence_files` (required, 1+): paths to tab-delimited sequence files, one per library
  - `--primer_file` (optional): path to user-defined primer file; defaults to `/data/pooled_primers.txt`

### 3.2 Additions to `pyguide/guide.py`

**New function: `read_gene_list_pooled_seq(file: str) -> Tuple[List, List, List, List, List]`**
- Reads the 5-column collated file: `name\tseq\tleft_primer\tright_primer\tlib_num`
- Returns `(names, seqs, left_primers, right_primers, lib_nums)`
- Mirrors the existing `read_gene_list_pooled` signature

**`order_guides()` — new `pooled-seq` branch:**
- Reads the collated file via `read_gene_list_pooled_seq`
- Builds `collated_df` directly with columns: `name`, `seq`, `gene` (set to `name`), `left_primer`, `right_primer`, `lib_num`
- Skips all database lookup
- Calls existing `write_pooled_txt(collated_df, ...)` unchanged
- Writes a simplified log file: date, user name, primer assignments per library (no gene alias / missing-gene sections)

**`main()` — additions:**
- Add `pooled-seq` to the `assert order_format in (...)` validity check
- Add `pooled-seq` to the primer-required guard (same as `pooled` and `batch-retest`)

### 3.3 `pyproject.toml`

Add one entry point:
```toml
pyguide-collate-seq = "pyguide.collate_seq:main"
```

---

## 4. File Formats

### 4.1 Input sequence file (user-created)

Tab-delimited, no header, one guide per line:

```
MY_GUIDE_1	ACGTACGTACGTACGTACGT
MY_GUIDE_2	TGCATGCATGCATGCATGCA
```

### 4.2 Collated output file (5-column, written by pyguide-collate-seq)

Tab-delimited, no header:

```
MY_GUIDE_1	ACGTACGTACGTACGTACGT	LEFT_PRIMER_SEQ	RIGHT_PRIMER_SEQ	0
MY_GUIDE_2	TGCATGCATGCATGCATGCA	LEFT_PRIMER_SEQ	RIGHT_PRIMER_SEQ	0
```

Columns: `name`, `seq`, `left_primer`, `right_primer`, `lib_num`

### 4.3 Validation error output (stderr)

```
Error: Invalid sequences in my_guides.txt:
  Row 3 "SHORT_GUIDE": length 4, expected 20
  Row 7 "BAD_CHARS":   "ACGTACGNACGTACGTACGT" contains invalid characters: N
Fix the above rows and try again.
```

---

## 5. Downstream pipeline (unchanged)

`write_pooled_txt` in `guide.py` already reads `collated_df['seq']` and wraps each spacer with:
- Left: `left_adapter + CCACCTTGTTG` (BstXI restriction site)
- Right: `GTTTAAGAGCTAAGCTGG + right_adapter` (Bpi1102I restriction site)
- Left and right PCR primers

No changes to `write_pooled_txt` or any other downstream function.

---

## 6. Testing

**Fixture file:** `testing/example/seq_guide_list.txt`
- 3–5 rows of `name\tsequence` using known-good STAT3 spacer sequences from the Horlbeck library

**New test file: `testing/test_seq.py`**

Unit tests:
- `test_validate_sequences_valid` — clean 20nt ACGT input passes without error
- `test_validate_sequences_wrong_length` — short/long sequence → raises with row number and name
- `test_validate_sequences_invalid_chars` — non-ACGT character → raises listing offending chars
- `test_validate_sequences_duplicate_names` — duplicate names → raises listing duplicates

End-to-end test:
- `test_collate_seq_and_order` — runs `collate_seq.main()` with fixture file, verifies collated file created; runs `guide.order_guides(..., order_format="pooled-seq", ...)`, verifies output `.txt` created; cleans up all generated files in `finally` block

---

## 7. Out of Scope

- CRISPRa vs CRISPRi distinction for raw sequences (the `--ai` flag is accepted for log consistency but has no effect on sequence processing — restriction sites and primers are identical for both)
- NTC generation for `pooled-seq` (NTCs are Horlbeck library sequences; not applicable to custom sequences)
- Single-guide or arrayed formats for raw sequences
- Sequence deduplication across libraries

---

## 8. Success Criteria

1. `pyguide-collate-seq --sequence_files my_guides.txt` produces a 5-column collated file
2. `pyguide-order --order_format pooled-seq` produces an Agilent-format `.txt` file identical in structure to existing `pooled` output
3. Invalid sequences produce clear error messages and a non-zero exit code
4. All existing tests continue to pass
5. New `test_seq.py` tests pass
