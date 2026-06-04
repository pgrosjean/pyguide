# Pooled-Seq Feature Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `pyguide-collate-seq` command and `pooled-seq` order format so users can supply raw 20nt spacer sequences directly (instead of gene names or Horlbeck guide IDs) and produce a pooled Agilent oligo order file.

**Architecture:** New `pyguide/collate_seq.py` handles validation and collation (reusing `get_primer_list` from `pool.py`). Two additions to `guide.py`: `read_gene_list_pooled_seq` reader and an early-return `pooled-seq` branch in `order_guides` that builds `collated_df` from `primer_df` directly (no database lookup) and calls the existing `write_pooled_txt` unchanged. A new `write_pooled_seq_log_file` writes a simplified log (primer assignments only).

**Tech Stack:** Python 3.12, pandas, uv, pytest, existing `pool.get_primer_list`, existing `guide.write_pooled_txt`

---

## File Map

| File | Action | What changes |
|------|--------|--------------|
| `pyguide/collate_seq.py` | **Create** | `validate_sequence_file`, `generate_pooled_seq_list`, `main` |
| `pyguide/guide.py` | **Modify** | Add `read_gene_list_pooled_seq`, `write_pooled_seq_log_file`, `pooled-seq` branch in `order_guides` and `main` |
| `pyproject.toml` | **Modify** | Add `pyguide-collate-seq` entry point |
| `testing/test_seq.py` | **Create** | Unit + end-to-end tests for the new feature |
| `testing/example/seq_guide_list.txt` | **Create** | Fixture: 3 STAT3 spacer sequences |

---

## Task 1: validate_sequence_file

**Files:**
- Create: `pyguide/collate_seq.py`
- Create: `testing/test_seq.py`

This function reads a tab-delimited sequence file (name, 20nt spacer), collects ALL validation errors before raising, and returns a clean DataFrame.

- [ ] **Step 1: Write the failing tests**

Create `testing/test_seq.py`:

```python
import pytest
import os
import glob
import pandas as pd
from pyguide import collate_seq, guide


def test_validate_sequences_valid(tmp_path):
    f = tmp_path / "guides.txt"
    f.write_text("MY_GUIDE_1\tACGTACGTACGTACGTACGT\nMY_GUIDE_2\tTGCATGCATGCATGCATGCA\n")
    df = collate_seq.validate_sequence_file(str(f))
    assert len(df) == 2
    assert list(df['name']) == ['MY_GUIDE_1', 'MY_GUIDE_2']
    assert list(df['seq']) == ['ACGTACGTACGTACGTACGT', 'TGCATGCATGCATGCATGCA']


def test_validate_sequences_wrong_length(tmp_path):
    f = tmp_path / "guides.txt"
    f.write_text("SHORT_GUIDE\tACGT\n")
    with pytest.raises(SystemExit):
        collate_seq.validate_sequence_file(str(f))


def test_validate_sequences_invalid_chars(tmp_path):
    f = tmp_path / "guides.txt"
    f.write_text("BAD_CHARS\tACGTACGNACGTACGTACGT\n")
    with pytest.raises(SystemExit):
        collate_seq.validate_sequence_file(str(f))


def test_validate_sequences_duplicate_names(tmp_path):
    f = tmp_path / "guides.txt"
    f.write_text("GUIDE_1\tACGTACGTACGTACGTACGT\nGUIDE_1\tTGCATGCATGCATGCATGCA\n")
    with pytest.raises(SystemExit):
        collate_seq.validate_sequence_file(str(f))
```

- [ ] **Step 2: Run tests to confirm they fail**

```bash
uv run pytest testing/test_seq.py -v
```

Expected: 4 errors — `ModuleNotFoundError: No module named 'pyguide.collate_seq'`

- [ ] **Step 3: Create `pyguide/collate_seq.py` with validate_sequence_file**

```python
import sys
import os
import pandas as pd
import numpy as np
from datetime import datetime
from typing import List, Tuple
from argparse import ArgumentParser
from pyguide.guide import get_unique_filename
from pyguide.pool import get_primer_list


def validate_sequence_file(file: str) -> pd.DataFrame:
    """
    Reads a tab-delimited sequence file (name, 20nt spacer) and validates every row.

    Collects ALL errors before exiting so the user sees every problem at once.
    Normalizes sequences to uppercase on success.

    Parameters
    ----------
    file : str
        Path to tab-delimited file with columns: name, sequence (no header).

    Returns
    -------
    pd.DataFrame with columns 'name' and 'seq' (uppercase).
    """
    df = pd.read_csv(file, sep="\t", header=None, names=['name', 'seq'])
    if len(df) == 0:
        print(f"Error: {file} is empty.", file=sys.stderr)
        sys.exit(1)
    df['seq'] = df['seq'].str.upper()
    errors = []
    seen_names: dict = {}
    for idx, row in df.iterrows():
        name, seq = str(row['name']), str(row['seq'])
        row_num = idx + 1
        if name in seen_names:
            errors.append(
                f"  Row {row_num} \"{name}\": duplicate name (first seen at row {seen_names[name]})"
            )
        else:
            seen_names[name] = row_num
        if len(seq) != 20:
            errors.append(f"  Row {row_num} \"{name}\": length {len(seq)}, expected 20")
        invalid_chars = set(seq) - {'A', 'C', 'G', 'T'}
        if invalid_chars:
            errors.append(
                f"  Row {row_num} \"{name}\": contains invalid characters: {''.join(sorted(invalid_chars))}"
            )
    if errors:
        print(f"Error: Invalid sequences in {file}:", file=sys.stderr)
        for error in errors:
            print(error, file=sys.stderr)
        print("Fix the above rows and try again.", file=sys.stderr)
        sys.exit(1)
    return df
```

- [ ] **Step 4: Run tests to confirm they pass**

```bash
uv run pytest testing/test_seq.py::test_validate_sequences_valid \
             testing/test_seq.py::test_validate_sequences_wrong_length \
             testing/test_seq.py::test_validate_sequences_invalid_chars \
             testing/test_seq.py::test_validate_sequences_duplicate_names -v
```

Expected: 4 PASS

- [ ] **Step 5: Run full suite to confirm no regressions**

```bash
uv run pytest -v
```

Expected: all tests pass.

- [ ] **Step 6: Commit**

```bash
git add pyguide/collate_seq.py testing/test_seq.py
git commit -m "$(cat <<'EOF'
feat: add validate_sequence_file for collate-seq pipeline

Validates tab-delimited name+sequence files before pooled ordering.
Collects all errors (wrong length, invalid chars, duplicate names)
before exiting so users see every problem at once.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: generate_pooled_seq_list + main() + entry point

**Files:**
- Modify: `pyguide/collate_seq.py`
- Modify: `pyproject.toml`

Adds the collation logic and CLI entry point. Reuses `get_primer_list` from `pool.py`.

- [ ] **Step 1: Add generate_pooled_seq_list and main() to collate_seq.py**

Append to `pyguide/collate_seq.py` (after `validate_sequence_file`):

```python
def generate_pooled_seq_list(sequence_files: List[str],
                             primer_list: List[Tuple[str, str]]) -> str:
    """
    Reads validated sequence files, assigns primers, and writes a 5-column
    collated file: name<TAB>seq<TAB>left_primer<TAB>right_primer<TAB>lib_num.

    Parameters
    ----------
    sequence_files : List[str]
        Paths to tab-delimited sequence files, one per library.
    primer_list : List[Tuple[str, str]]
        Primer pairs, one per library (same length as sequence_files).

    Returns
    -------
    str
        Path to the written collated file.
    """
    assert len(sequence_files) == len(primer_list), \
        "sequence_files and primer_list must have the same length."
    base_dir = os.path.dirname(os.path.abspath(sequence_files[0]))
    now = datetime.now()
    date = now.strftime("%y_%m_%d")
    file_name = f"collated_seq_wishlist_{date}.txt"
    file_name = get_unique_filename(base_dir, file_name)
    out_path = os.path.join(base_dir, file_name)
    with open(out_path, 'w') as out_file:
        for lib_num, (seq_file, primer_tuple) in enumerate(zip(sequence_files, primer_list)):
            df = validate_sequence_file(seq_file)
            left_primer, right_primer = primer_tuple
            for _, row in df.iterrows():
                out_file.write(f"{row['name']}\t{row['seq']}\t{left_primer}\t{right_primer}\t{lib_num}\n")
    print(f"Collated sequence wishlist generated: {out_path}")
    return out_path


def main(raw_args=None):
    parser = ArgumentParser()
    parser.add_argument("--sequence_files",
                        type=str,
                        nargs="+",
                        required=True,
                        help="Paths to tab-delimited sequence files (name<TAB>20nt spacer), one per library.")
    parser.add_argument("--primer_file",
                        type=str,
                        default=None,
                        help="Path to user-defined primer file. Defaults to /data/pooled_primers.txt.")
    args = parser.parse_args(raw_args)

    file_path = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    default_primer_file = os.path.join(file_path, '..', 'data', 'pooled_primers.txt')
    primer_file = args.primer_file if args.primer_file is not None else default_primer_file

    num_pairs = len(args.sequence_files)
    primer_list = get_primer_list(primer_file, num_pairs)
    generate_pooled_seq_list(args.sequence_files, primer_list)


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Add entry point to pyproject.toml**

In `pyproject.toml`, find the `[project.scripts]` section:

```toml
[project.scripts]
pyguide-order = "pyguide.guide:main"
pyguide-collate = "pyguide.pool:main"
pyguide-batch-retest = "pyguide.batch_retest:main"
pyguide-check-seq = "pyguide.check_seq:main"
```

Replace with:

```toml
[project.scripts]
pyguide-order = "pyguide.guide:main"
pyguide-collate = "pyguide.pool:main"
pyguide-collate-seq = "pyguide.collate_seq:main"
pyguide-batch-retest = "pyguide.batch_retest:main"
pyguide-check-seq = "pyguide.check_seq:main"
```

- [ ] **Step 3: Re-sync so the new entry point is installed**

```bash
uv sync --extra dev
```

Expected: no errors.

- [ ] **Step 4: Verify the entry point is callable**

```bash
uv run pyguide-collate-seq --help
```

Expected: prints usage with `--sequence_files` and `--primer_file` options.

- [ ] **Step 5: Run full suite**

```bash
uv run pytest -v
```

Expected: all tests pass.

- [ ] **Step 6: Commit**

```bash
git add pyguide/collate_seq.py pyproject.toml uv.lock
git commit -m "$(cat <<'EOF'
feat: add pyguide-collate-seq command

Reads tab-delimited name+sequence files, assigns primers from the
default pooled_primers.txt pool (or user-supplied file), and writes
a 5-column collated_seq_wishlist_<date>.txt for use with
pyguide-order --order_format pooled-seq.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: read_gene_list_pooled_seq in guide.py

**Files:**
- Modify: `pyguide/guide.py` (after `read_gene_list_pooled`, around line 100)
- Modify: `testing/test_seq.py`

Adds the reader for the 5-column collated file. Mirrors `read_gene_list_pooled`.

- [ ] **Step 1: Write the failing test**

Add to `testing/test_seq.py`:

```python
def test_read_gene_list_pooled_seq(tmp_path):
    f = tmp_path / "collated.txt"
    f.write_text(
        "GUIDE_A\tACGTACGTACGTACGTACGT\tLEFT1\tRIGHT1\t0\n"
        "GUIDE_B\tTGCATGCATGCATGCATGCA\tLEFT1\tRIGHT1\t0\n"
    )
    names, seqs, lefts, rights, lib_nums = guide.read_gene_list_pooled_seq(str(f))
    assert names == ['GUIDE_A', 'GUIDE_B']
    assert seqs == ['ACGTACGTACGTACGTACGT', 'TGCATGCATGCATGCATGCA']
    assert lefts == ['LEFT1', 'LEFT1']
    assert rights == ['RIGHT1', 'RIGHT1']
    assert lib_nums == [0, 0]
```

- [ ] **Step 2: Run test to confirm it fails**

```bash
uv run pytest testing/test_seq.py::test_read_gene_list_pooled_seq -v
```

Expected: FAIL — `AttributeError: module 'pyguide.guide' has no attribute 'read_gene_list_pooled_seq'`

- [ ] **Step 3: Add read_gene_list_pooled_seq to guide.py**

In `pyguide/guide.py`, immediately after `read_gene_list_pooled` (after line 99), insert:

```python
def read_gene_list_pooled_seq(file: str) -> Tuple[List[str], List[str], List[str], List[str], List]:
    """
    Reads a 5-column collated sequence wishlist file produced by pyguide-collate-seq.

    Parameters
    ----------
    file : str
        Path to the collated file (name, seq, left_primer, right_primer, lib_num).

    Returns
    -------
    names, seqs, left_primers, right_primers, lib_nums
    """
    df = pd.read_csv(file, sep="\t", header=None)
    assert len(df.columns) == 5, \
        f"Expected 5-column file from pyguide-collate-seq. Got {len(df.columns)} columns. {df}"
    names = list(df[0].values)
    seqs = list(df[1].values)
    left_primers = list(df[2].values)
    right_primers = list(df[3].values)
    lib_nums = list(df[4].values)
    return names, seqs, left_primers, right_primers, lib_nums
```

- [ ] **Step 4: Run test to confirm it passes**

```bash
uv run pytest testing/test_seq.py::test_read_gene_list_pooled_seq -v
```

Expected: PASS

- [ ] **Step 5: Run full suite**

```bash
uv run pytest -v
```

Expected: all tests pass.

- [ ] **Step 6: Commit**

```bash
git add pyguide/guide.py testing/test_seq.py
git commit -m "$(cat <<'EOF'
feat: add read_gene_list_pooled_seq reader to guide.py

Reads the 5-column collated_seq_wishlist file produced by
pyguide-collate-seq. Mirrors read_gene_list_pooled but returns
an extra sequence column.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: write_pooled_seq_log_file in guide.py

**Files:**
- Modify: `pyguide/guide.py` (after `write_batch_retest_log_file`, around line 993)
- Modify: `testing/test_seq.py`

Adds the log writer for pooled-seq ordering. Mirrors `write_batch_retest_log_file`.

- [ ] **Step 1: Write the failing test**

Add to `testing/test_seq.py`:

```python
def test_write_pooled_seq_log_file(tmp_path):
    primer_df = pd.DataFrame({
        'guide_id': ['G1', 'G2'],
        'seq': ['ACGTACGTACGTACGTACGT', 'TGCATGCATGCATGCATGCA'],
        'left_primers': ['LEFTSEQ', 'LEFTSEQ'],
        'right_primers': ['RIGHTSEQ', 'RIGHTSEQ'],
        'lib_num': [0, 0],
    })
    guide.write_pooled_seq_log_file("TestUser", str(tmp_path), primer_df)
    log_files = list(tmp_path.glob("log_file_pooled_seq_TestUser_*.txt"))
    assert len(log_files) == 1
    content = log_files[0].read_text()
    assert "TestUser" in content
    assert "LEFTSEQ" in content
    assert "RIGHTSEQ" in content
```

- [ ] **Step 2: Run test to confirm it fails**

```bash
uv run pytest testing/test_seq.py::test_write_pooled_seq_log_file -v
```

Expected: FAIL — `AttributeError: module 'pyguide.guide' has no attribute 'write_pooled_seq_log_file'`

- [ ] **Step 3: Add write_pooled_seq_log_file to guide.py**

In `pyguide/guide.py`, immediately after `write_batch_retest_log_file` (after line 992), insert:

```python
def write_pooled_seq_log_file(name: str, base_dir: str, primer_df: pd.DataFrame):
    """
    Writes a log file for pooled-seq ordering. Records date, user name,
    and primer assignments per library. No gene-alias or missing-gene
    sections (sequences are supplied directly).

    Parameters
    ----------
    name : str
        Name of user.
    base_dir : str
        Directory to save the log file.
    primer_df : pd.DataFrame
        DataFrame with columns guide_id, seq, left_primers, right_primers, lib_num.
    """
    date = get_current_date()
    filename = get_unique_filename(base_dir, f"log_file_pooled_seq_{name}_{date}.txt")
    file = os.path.join(base_dir, filename)
    unique_primers = primer_df[['left_primers', 'right_primers', 'lib_num']].drop_duplicates()
    with open(file, 'w') as f:
        f.write(f"This pooled-seq library was generated by {name} on {date} \n \n")
        f.write(f"\nPrimer Information\n------------------\nLibrary Number\tPrimer 1\tPrimer 2\n")
        for _, row in unique_primers.iterrows():
            f.write(f"{row['lib_num']}\t{row['left_primers']}\t{row['right_primers']}\n")
```

- [ ] **Step 4: Run test to confirm it passes**

```bash
uv run pytest testing/test_seq.py::test_write_pooled_seq_log_file -v
```

Expected: PASS

- [ ] **Step 5: Run full suite**

```bash
uv run pytest -v
```

Expected: all tests pass.

- [ ] **Step 6: Commit**

```bash
git add pyguide/guide.py testing/test_seq.py
git commit -m "$(cat <<'EOF'
feat: add write_pooled_seq_log_file to guide.py

Writes a simplified pooled-seq log (date, user, primer assignments).
No gene-alias or missing-gene sections since sequences are supplied
directly rather than looked up from the Horlbeck database.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: pooled-seq branch in order_guides() and main()

**Files:**
- Modify: `pyguide/guide.py`

Wires the pooled-seq pathway end-to-end inside the existing CLI. The early-return pattern keeps this completely isolated from all existing code paths.

- [ ] **Step 1: Add pooled-seq to order_guides()**

In `pyguide/guide.py`, find these two lines near the top of `order_guides` (around line 1135):

```python
    possible_order_formats = ['arrayed', 'single', 'pooled', 'batch-retest']
    assert order_format in possible_order_formats, f"{order_format} not in {possible_order_formats}"
    assert guides_per_gene <= 10, "No more than 10 guides per gene allowed."
```

Replace with:

```python
    possible_order_formats = ['arrayed', 'single', 'pooled', 'batch-retest', 'pooled-seq']
    assert order_format in possible_order_formats, f"{order_format} not in {possible_order_formats}"
    assert guides_per_gene <= 10, "No more than 10 guides per gene allowed."

    if order_format == "pooled-seq":
        assert primer_df is not None, \
            "You must run pyguide-collate-seq and provide the collated file before pooled-seq ordering."
        collated_df = pd.DataFrame({
            'name': primer_df['guide_id'].values,
            'seq': primer_df['seq'].values,
            'left_primers': primer_df['left_primers'].values,
            'right_primers': primer_df['right_primers'].values,
            'lib_num': primer_df['lib_num'].values,
            'gene': primer_df['guide_id'].values,
        })
        write_pooled_txt(collated_df, name, base_dir)
        write_pooled_seq_log_file(name, base_dir, primer_df)
        return
```

- [ ] **Step 2: Add pooled-seq to main()**

In `pyguide/guide.py` `main()`, find the `if args.order_format.lower() == "pooled":` block (around line 1340) and add a new `elif` for `pooled-seq` immediately after the `batch-retest` block:

Find:
```python
    else:
        guide_list, gene_list = read_wishlist(file=file)
        primer_df = None
```

Replace with:
```python
    elif args.order_format.lower() == "pooled-seq":
        name_list, seq_list, left_primers, right_primers, lib_num = read_gene_list_pooled_seq(file=file)
        guide_list = []
        gene_list = []
        primer_df = pd.DataFrame({
            'guide_id': name_list,
            'seq': seq_list,
            'left_primers': left_primers,
            'right_primers': right_primers,
            'lib_num': lib_num,
        })
    else:
        guide_list, gene_list = read_wishlist(file=file)
        primer_df = None
```

- [ ] **Step 3: Update the order_format assertion and primer guard in main()**

Find (around line 1368):
```python
    assert args.order_format in ["single", "pooled", "arrayed", "batch-retest"], "Only single, pooled, batch-retest, or arrayed order formats."
    if args.order_format == "pooled" or args.order_format == "batch-retest":
        assert primer_df is not None, "You must run pyguide-collate (or pyguide-batch-retest) and generate primers before pooled ordering."
```

Replace with:
```python
    assert args.order_format in ["single", "pooled", "pooled-seq", "arrayed", "batch-retest"], \
        "Only single, pooled, pooled-seq, arrayed, or batch-retest order formats."
    if args.order_format in ("pooled", "batch-retest", "pooled-seq"):
        assert primer_df is not None, \
            "You must run pyguide-collate (or pyguide-collate-seq or pyguide-batch-retest) before pooled ordering."
```

- [ ] **Step 4: Run full suite**

```bash
uv run pytest -v
```

Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add pyguide/guide.py
git commit -m "$(cat <<'EOF'
feat: add pooled-seq order format to order_guides and main

Early-return branch in order_guides builds collated_df directly
from primer_df (which carries both sequences and primer assignments),
bypassing the Horlbeck database lookup entirely. write_pooled_txt
is called unchanged so restriction sites and adapters are added
identically to the standard pooled format.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Fixture file and end-to-end test

**Files:**
- Create: `testing/example/seq_guide_list.txt`
- Modify: `testing/test_seq.py`

Adds the fixture and full pipeline test.

- [ ] **Step 1: Create the fixture file**

Create `testing/example/seq_guide_list.txt` with exactly this content (3 real STAT3 spacer sequences from the Horlbeck CRISPRi library):

```
STAT3_guide_1	GCCAGGGAGCAGGAAATCGG
STAT3_guide_2	GGGATCCTGGACAGGCACCC
STAT3_guide_3	GAGGGCCTCTCCGAGCCGAG
```

- [ ] **Step 2: Write the end-to-end test**

Add to `testing/test_seq.py`:

```python
def test_collate_seq_and_order():
    file_path_1 = os.path.dirname(os.path.abspath(__file__))
    seq_file = os.path.join(file_path_1, "example", "seq_guide_list.txt")
    example_dir = os.path.join(file_path_1, "example")

    try:
        # Step 1: collate
        collate_seq.main(["--sequence_files", seq_file])
        collated = glob.glob(os.path.join(example_dir, "collated_seq_wishlist_*.txt"))
        assert len(collated) >= 1, "collate_seq produced no collated file"
        collated_file = sorted(collated, key=os.path.getmtime)[-1]

        # Step 2: read collated file and build primer_df
        names, seqs, left_primers, right_primers, lib_nums = guide.read_gene_list_pooled_seq(collated_file)
        primer_df = pd.DataFrame({
            'guide_id': names,
            'seq': seqs,
            'left_primers': left_primers,
            'right_primers': right_primers,
            'lib_num': lib_nums,
        })

        # Step 3: order
        guide.order_guides(
            guide_ids=[],
            gene_names=[],
            name="Test",
            ai_status="i",
            guides_per_gene=5,
            order_format="pooled-seq",
            base_dir=example_dir,
            check_db=False,
            organism="human",
            primer_df=primer_df,
        )

        # Step 4: verify output file exists and has content
        order_files = glob.glob(os.path.join(example_dir, "order_pooled_Test_*.txt"))
        assert len(order_files) >= 1, "order_guides(pooled-seq) produced no output TXT"
        with open(sorted(order_files, key=os.path.getmtime)[-1]) as fh:
            lines = [l for l in fh.readlines() if l.strip()]
        assert len(lines) >= 6, f"Expected at least 6 lines (2 per guide × 3 guides), got {len(lines)}"

    finally:
        for f in (
            glob.glob(os.path.join(example_dir, "order_pooled_Test_*.txt"))
            + glob.glob(os.path.join(example_dir, "order_pooled_Test_*_info.csv"))
            + glob.glob(os.path.join(example_dir, "log_file_pooled_seq_Test_*.txt"))
            + glob.glob(os.path.join(example_dir, "collated_seq_wishlist_*.txt"))
        ):
            if os.path.exists(f):
                os.remove(f)
```

- [ ] **Step 3: Run the end-to-end test**

```bash
uv run pytest testing/test_seq.py::test_collate_seq_and_order -v
```

Expected: PASS

- [ ] **Step 4: Run full suite**

```bash
uv run pytest -v
```

Expected: all tests pass.

- [ ] **Step 5: Verify the CLI works end-to-end**

```bash
uv run pyguide-collate-seq --help
```

Expected: prints usage without errors.

- [ ] **Step 6: Commit**

```bash
git add testing/test_seq.py testing/example/seq_guide_list.txt
git commit -m "$(cat <<'EOF'
test: add end-to-end tests for pooled-seq pipeline

Unit tests for validate_sequence_file (valid input, wrong length,
invalid chars, duplicate names), read_gene_list_pooled_seq,
write_pooled_seq_log_file, and a full pipeline test that runs
collate_seq → order_guides(pooled-seq) and verifies output.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Final Verification

- [ ] **Run full test suite:**

```bash
uv run pytest -v
```

Expected: all tests pass (original 20 + new tests from test_seq.py).

- [ ] **Verify all entry points:**

```bash
uv run pyguide-collate-seq --help
uv run pyguide-order --help
```

Expected: both print help text. `pyguide-order --help` shows `pooled-seq` is an accepted `--order_format` value (via the assert message).

- [ ] **Check git log:**

```bash
git log --oneline parker/update-infra ^master
```

Expected: commits covering spec, README, all 11 bug fixes, test improvements, and pooled-seq feature (Tasks 1–6 above plus earlier infra work).
