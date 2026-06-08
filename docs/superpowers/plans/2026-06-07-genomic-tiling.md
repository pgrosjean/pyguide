# Genomic Tiling Library Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `pyguide-tiling`, a CLI command that takes chromosomal coordinates, runs GuideScan2 to find all NGG-PAM SpCas9 guides, filters them, and outputs a tab-delimited sequence file compatible with the `pyguide-collate-seq` → `pyguide-order --order_format pooled-seq` pipeline.

**Architecture:** User provides coordinates (e.g. `chr1:12345-12395`); the tool fetches the genomic sequence via UCSC REST API, scans for all NGG PAM sites on both strands, runs `guidescan enumerate` to get per-guide specificity scores, applies filters (TTTT motif, restriction sites, specificity threshold), then greedily selects non-redundant guides by pairwise Hamming distance (prioritizing highest specificity), and optionally caps guides per region with even positional spread. Optionally writes a BigWig coverage track.

**Tech Stack:** Python 3.12, pandas, numpy, subprocess (guidescan CLI), urllib.request (UCSC REST API), pyBigWig (optional), hatchling

---

## File Structure

**New files:**
- `pyguide/tiling.py` — all pipeline logic + `main()` entry point
- `pyguide/tiling_scripts/__init__.py` — empty, makes it a package
- `pyguide/tiling_scripts/decode_database.py` — copied from guidescan-cli (MIT)
- `pyguide/tiling_scripts/cfd/mismatch_score.pkl` — copied from guidescan-cli (MIT)
- `pyguide/tiling_scripts/cfd/pam_scores.pkl` — copied from guidescan-cli (MIT)
- `pyguide/data/hg38.chrom.sizes` — hg38 chromosome sizes for BigWig
- `scripts/download_guidescan_data.sh` — downloads GuideScan2 index (~2.2 GB)
- `testing/test_tiling.py` — all unit tests

**Modified files:**
- `pyproject.toml` — add `pyguide-tiling` entry point + `tiling = ["pyBigWig"]` optional dep
- `README.md` — GuideScan2 setup, tiling usage, BigWig output

---

## Key Design Notes (read before implementing)

**GuideScan2 CLI interface:**
`guidescan enumerate` takes a kmers file (`-f`/`--kmers-file`, one 20nt sequence per line), an index (`--index`), and writes CSV output to a file (`-o`). The CSV header is:
```
id,sequence,match_chrm,match_position,match_strand,match_distance,specificity
```
Each row is one match (on-target + off-targets). We group by `id` and take the first `specificity` value per guide (it's the same for all rows of the same guide — it's an aggregate property of the guide, not per-match).

**Coordinate convention:** User input is 1-based closed (`chr1:12345-12395`). UCSC REST API takes 0-based half-open (pass `start-1, end`). Guide positions in our output use 1-based.

**Sequence fetching:** We use the UCSC REST API (`https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=X;start=Y;end=Z`) to get the region sequence. No local FASTA needed.

**GuideScan2 not on PyPI.** At runtime, check `shutil.which("guidescan")`. If None, raise `RuntimeError` with install instructions. Tests mock this check.

---

## Task 1: Branch + pyproject.toml

**Files:**
- Modify: `pyproject.toml`

- [ ] **Step 1: Create the genomic_tiling branch**

```bash
git checkout -b genomic_tiling
```

- [ ] **Step 2: Write the failing test (verify pyproject changes take effect)**

In `testing/test_tiling.py`:
```python
import pytest

def test_import_tiling():
    from pyguide import tiling
    assert hasattr(tiling, 'parse_coordinates')
```

- [ ] **Step 3: Run test to confirm it fails**

```bash
cd /path/to/pyguide && uv run pytest testing/test_tiling.py::test_import_tiling -v
```
Expected: FAIL with `ModuleNotFoundError` or `ImportError`

- [ ] **Step 4: Update pyproject.toml**

Replace the `[project.optional-dependencies]` and `[project.scripts]` sections:
```toml
[project.optional-dependencies]
dev = ["pytest>=7.0"]
app = ["streamlit"]
tiling = ["pyBigWig"]

[project.scripts]
pyguide-order = "pyguide.guide:main"
pyguide-collate = "pyguide.pool:main"
pyguide-collate-seq = "pyguide.collate_seq:main"
pyguide-batch-retest = "pyguide.batch_retest:main"
pyguide-check-seq = "pyguide.check_seq:main"
pyguide-tiling = "pyguide.tiling:main"
```

- [ ] **Step 5: Create `pyguide/tiling.py` with a stub**

```python
def parse_coordinates(coord_str: str) -> tuple[str, int, int]:
    raise NotImplementedError


def main(raw_args=None):
    pass
```

- [ ] **Step 6: Run test to confirm it passes**

```bash
uv run pytest testing/test_tiling.py::test_import_tiling -v
```
Expected: PASS

- [ ] **Step 7: Commit**

```bash
git add pyproject.toml pyguide/tiling.py testing/test_tiling.py
git commit -m "feat: scaffold pyguide-tiling entry point and tiling optional dep"
```

---

## Task 2: Bundle GuideScan2 CFD scripts and hg38.chrom.sizes

**Files:**
- Create: `pyguide/tiling_scripts/__init__.py`
- Create: `pyguide/tiling_scripts/decode_database.py`
- Create: `pyguide/tiling_scripts/cfd/mismatch_score.pkl`
- Create: `pyguide/tiling_scripts/cfd/pam_scores.pkl`
- Create: `pyguide/data/hg38.chrom.sizes`

- [ ] **Step 1: Create directories**

```bash
mkdir -p pyguide/tiling_scripts/cfd
mkdir -p pyguide/data
touch pyguide/tiling_scripts/__init__.py
```

- [ ] **Step 2: Download decode_database.py from guidescan-cli (MIT license)**

```bash
curl -o pyguide/tiling_scripts/decode_database.py \
  https://raw.githubusercontent.com/pritykinlab/guidescan-cli/master/scripts/decode_database.py
```

- [ ] **Step 3: Download CFD pickle files**

```bash
curl -o pyguide/tiling_scripts/cfd/mismatch_score.pkl \
  https://raw.githubusercontent.com/pritykinlab/guidescan-cli/master/scripts/cfd/mismatch_score.pkl

curl -o pyguide/tiling_scripts/cfd/pam_scores.pkl \
  https://raw.githubusercontent.com/pritykinlab/guidescan-cli/master/scripts/cfd/pam_scores.pkl
```

- [ ] **Step 4: Download hg38.chrom.sizes from UCSC**

```bash
curl -o pyguide/data/hg38.chrom.sizes \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
```

Verify the file exists and has content:
```bash
head -5 pyguide/data/hg38.chrom.sizes
```
Expected output (first few lines):
```
chr1	248956422
chr2	242193529
chr3	198295559
```

- [ ] **Step 5: Verify CFD pkl files are valid Python pickles**

```bash
uv run python -c "
import pickle
mm = pickle.load(open('pyguide/tiling_scripts/cfd/mismatch_score.pkl', 'rb'))
pam = pickle.load(open('pyguide/tiling_scripts/cfd/pam_scores.pkl', 'rb'))
print('mismatch keys sample:', list(mm.keys())[:3])
print('pam keys:', list(pam.keys()))
"
```
Expected: prints some keys without error.

- [ ] **Step 6: Commit**

```bash
git add pyguide/tiling_scripts/ pyguide/data/hg38.chrom.sizes
git commit -m "feat: bundle guidescan-cli CFD scripts and hg38 chrom sizes (MIT)"
```

---

## Task 3: parse_coordinates, sequence utilities, find_ngg_guides

**Files:**
- Modify: `pyguide/tiling.py`
- Modify: `testing/test_tiling.py`

- [ ] **Step 1: Write the failing tests**

Add to `testing/test_tiling.py`:
```python
import pandas as pd
from pyguide.tiling import (
    parse_coordinates,
    reverse_complement,
    fetch_sequence_ucsc,
    find_ngg_guides,
)


def test_parse_coordinates_valid():
    chrom, start, end = parse_coordinates("chr1:12345-12395")
    assert chrom == "chr1"
    assert start == 12345
    assert end == 12395


def test_parse_coordinates_invalid_format():
    with pytest.raises(ValueError, match="Expected format"):
        parse_coordinates("12345-12395")


def test_parse_coordinates_invalid_range():
    with pytest.raises(ValueError, match="start must be less than end"):
        parse_coordinates("chr1:12395-12345")


def test_reverse_complement():
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AAAA") == "TTTT"
    assert reverse_complement("GCTA") == "TAGC"


def test_find_ngg_guides_forward():
    # Sequence with a single NGG at position 20-22 (0-indexed)
    # spacer = first 20 nt, PAM = AGG
    seq = "ACGTACGTACGTACGTACGTAGG"  # 23nt: 20nt spacer + AGG
    guides = find_ngg_guides("chr1", 1, seq)
    fwd = [g for g in guides if g['match_strand'] == '+']
    assert len(fwd) >= 1
    assert fwd[0]['sequence'] == "ACGTACGTACGTACGTACGT"
    assert fwd[0]['match_strand'] == '+'
    assert fwd[0]['match_position'] == 1


def test_find_ngg_guides_reverse():
    # CCN on + strand means NGG guide on - strand
    # CC at positions 0-1 means spacer on - strand = revcomp(seq[3:23])
    seq = "CCT" + "ACGTACGTACGTACGTACGT"  # 23nt: CCN + 20nt
    guides = find_ngg_guides("chr1", 1, seq)
    rev = [g for g in guides if g['match_strand'] == '-']
    assert len(rev) >= 1
    expected_spacer = reverse_complement("ACGTACGTACGTACGTACGT")
    assert rev[0]['sequence'] == expected_spacer
    assert rev[0]['match_strand'] == '-'


def test_find_ngg_guides_excludes_n():
    # Spacer with N in it should be excluded
    seq = "ACGTACGTACGTACGNACGTAGG"
    guides = find_ngg_guides("chr1", 1, seq)
    assert not any('N' in g['sequence'] for g in guides)
```

- [ ] **Step 2: Run tests to confirm they fail**

```bash
uv run pytest testing/test_tiling.py::test_parse_coordinates_valid \
  testing/test_tiling.py::test_reverse_complement \
  testing/test_tiling.py::test_find_ngg_guides_forward -v
```
Expected: FAIL (NotImplementedError or ImportError)

- [ ] **Step 3: Implement utilities in pyguide/tiling.py**

Replace the stub with:
```python
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.request
import warnings
from argparse import ArgumentParser
from typing import Optional

import pandas as pd


_COMP = str.maketrans('ACGTN', 'TGCAN')


def reverse_complement(seq: str) -> str:
    return seq.upper().translate(_COMP)[::-1]


def parse_coordinates(coord_str: str) -> tuple[str, int, int]:
    match = re.fullmatch(r'(chr\w+):(\d+)-(\d+)', coord_str.strip())
    if not match:
        raise ValueError(
            f"Invalid coordinate '{coord_str}'. Expected format: chr<N>:start-end (e.g. chr1:12345-12395)"
        )
    chrom = match.group(1)
    start = int(match.group(2))
    end = int(match.group(3))
    if start >= end:
        raise ValueError(
            f"Invalid coordinate '{coord_str}': start must be less than end."
        )
    return chrom, start, end


def fetch_sequence_ucsc(chrom: str, start: int, end: int) -> str:
    """Fetch genomic sequence from UCSC REST API. start/end are 1-based closed."""
    # UCSC API uses 0-based, half-open coordinates
    url = (
        f"https://api.genome.ucsc.edu/getData/sequence"
        f"?genome=hg38;chrom={chrom};start={start - 1};end={end}"
    )
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            data = json.loads(resp.read())
    except Exception as e:
        raise RuntimeError(
            f"Failed to fetch sequence for {chrom}:{start}-{end} from UCSC: {e}"
        )
    return data['dna'].upper()


def find_ngg_guides(chrom: str, start: int, sequence: str) -> list[dict]:
    """
    Find all 20nt spacers with NGG PAM in sequence.

    start: 1-based genomic coordinate of sequence[0].
    Returns list of dicts: sequence, match_chrm, match_position (1-based), match_strand.
    """
    guides = []
    n = len(sequence)

    # Forward (+) strand: spacer at [i:i+20], PAM at [i+20:i+23] must be xGG
    for i in range(n - 22):
        if sequence[i + 21] == 'G' and sequence[i + 22] == 'G':
            spacer = sequence[i:i + 20]
            if 'N' not in spacer and set(spacer).issubset('ACGT'):
                guides.append({
                    'sequence': spacer,
                    'match_chrm': chrom,
                    'match_position': start + i,  # 1-based start of spacer
                    'match_strand': '+',
                })

    # Reverse (-) strand: CCx on + strand at [i:i+3], spacer = revcomp([i+3:i+23])
    for i in range(n - 22):
        if sequence[i] == 'C' and sequence[i + 1] == 'C':
            target = sequence[i + 3:i + 23]
            spacer = reverse_complement(target)
            if 'N' not in spacer and set(spacer).issubset('ACGT'):
                # 1-based 5' position of spacer on - strand = 3' end on + strand
                guides.append({
                    'sequence': spacer,
                    'match_chrm': chrom,
                    'match_position': start + i + 22,
                    'match_strand': '-',
                })

    return guides
```

- [ ] **Step 4: Run tests to confirm they pass**

```bash
uv run pytest testing/test_tiling.py::test_parse_coordinates_valid \
  testing/test_tiling.py::test_parse_coordinates_invalid_format \
  testing/test_tiling.py::test_parse_coordinates_invalid_range \
  testing/test_tiling.py::test_reverse_complement \
  testing/test_tiling.py::test_find_ngg_guides_forward \
  testing/test_tiling.py::test_find_ngg_guides_reverse \
  testing/test_tiling.py::test_find_ngg_guides_excludes_n -v
```
Expected: all PASS

- [ ] **Step 5: Commit**

```bash
git add pyguide/tiling.py testing/test_tiling.py
git commit -m "feat: add parse_coordinates, fetch_sequence_ucsc, find_ngg_guides"
```

---

## Task 4: run_guidescan — subprocess + specificity lookup

**Files:**
- Modify: `pyguide/tiling.py`
- Modify: `testing/test_tiling.py`

- [ ] **Step 1: Write the failing tests**

Add to `testing/test_tiling.py`:
```python
import shutil
from unittest.mock import patch, MagicMock
from pyguide.tiling import run_guidescan


def test_run_guidescan_missing_binary():
    guides = [
        {'sequence': 'ACGTACGTACGTACGTACGT', 'match_chrm': 'chr1',
         'match_position': 100, 'match_strand': '+'},
    ]
    with patch('shutil.which', return_value=None):
        with pytest.raises(RuntimeError, match="guidescan.*not found"):
            run_guidescan(guides, index_path="/fake/index")


def test_run_guidescan_returns_dataframe():
    """run_guidescan merges specificity from CSV onto guide list."""
    guides = [
        {'sequence': 'ACGTACGTACGTACGTACGT', 'match_chrm': 'chr1',
         'match_position': 100, 'match_strand': '+'},
    ]
    fake_csv = (
        "id,sequence,match_chrm,match_position,match_strand,match_distance,specificity\n"
        "guide_0,ACGTACGTACGTACGTACGT,chr1,99,+,0,0.85\n"
    )

    def fake_subprocess(cmd, **kwargs):
        # Write what guidescan would write to --output
        out_idx = cmd.index('--output')
        with open(cmd[out_idx + 1], 'w') as f:
            f.write(fake_csv)
        m = MagicMock()
        m.returncode = 0
        m.stderr = ''
        return m

    with patch('shutil.which', return_value='/usr/bin/guidescan'):
        with patch('subprocess.run', side_effect=fake_subprocess):
            df = run_guidescan(guides, index_path="/fake/index")
    assert 'specificity' in df.columns
    assert len(df) == 1
    assert df.iloc[0]['specificity'] == pytest.approx(0.85)
```

- [ ] **Step 2: Run tests to confirm they fail**

```bash
uv run pytest testing/test_tiling.py::test_run_guidescan_missing_binary \
  testing/test_tiling.py::test_run_guidescan_returns_dataframe -v
```
Expected: FAIL (NotImplementedError)

- [ ] **Step 3: Implement run_guidescan in pyguide/tiling.py**

Add after `find_ngg_guides`:
```python
def run_guidescan(guides: list[dict], index_path: str) -> pd.DataFrame:
    """
    Run guidescan enumerate on a list of guide dicts and return a DataFrame
    with all guide fields plus a 'specificity' column.

    guides: list of dicts with keys: sequence, match_chrm, match_position, match_strand
    index_path: path to the GuideScan2 hg38 index file
    """
    if shutil.which("guidescan") is None:
        raise RuntimeError(
            "guidescan not found on PATH. Install GuideScan2 (e.g. via Conda: "
            "`conda install -c bioconda guidescan2`) and ensure it is on your PATH. "
            "See README for full setup instructions."
        )

    if not guides:
        return pd.DataFrame(columns=[
            'sequence', 'match_chrm', 'match_position', 'match_strand', 'specificity'
        ])

    guides_df = pd.DataFrame(guides)

    with tempfile.TemporaryDirectory() as tmpdir:
        kmers_file = os.path.join(tmpdir, "kmers.txt")
        output_file = os.path.join(tmpdir, "output.csv")

        # Write one 20nt spacer per line
        with open(kmers_file, 'w') as f:
            for seq in guides_df['sequence']:
                f.write(seq + '\n')

        result = subprocess.run(
            [
                "guidescan", "enumerate",
                "--index", index_path,
                "--kmers-file", kmers_file,
                "--format", "csv",
                "--output", output_file,
            ],
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"guidescan enumerate failed:\n{result.stderr}"
            )

        if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
            # No guides found in index — return guides with NaN specificity
            guides_df['specificity'] = float('nan')
            return guides_df

        gs_df = pd.read_csv(output_file)

    # GuideScan2 may return multiple rows per guide (one per match).
    # Take the first specificity value per unique sequence (it's an aggregate property).
    if 'specificity' not in gs_df.columns:
        guides_df['specificity'] = float('nan')
        return guides_df

    spec_map = (
        gs_df.dropna(subset=['specificity'])
             .groupby('sequence')['specificity']
             .first()
             .to_dict()
    )

    guides_df['specificity'] = guides_df['sequence'].map(spec_map)
    return guides_df
```

- [ ] **Step 4: Run tests to confirm they pass**

```bash
uv run pytest testing/test_tiling.py::test_run_guidescan_missing_binary \
  testing/test_tiling.py::test_run_guidescan_returns_dataframe -v
```
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add pyguide/tiling.py testing/test_tiling.py
git commit -m "feat: add run_guidescan with subprocess and specificity merge"
```

---

## Task 5: apply_filters

**Files:**
- Modify: `pyguide/tiling.py`
- Modify: `testing/test_tiling.py`

- [ ] **Step 1: Write the failing tests**

Add to `testing/test_tiling.py`:
```python
from pyguide.tiling import apply_filters


def _make_guide_df(sequences, specificities=None):
    if specificities is None:
        specificities = [0.9] * len(sequences)
    return pd.DataFrame({
        'sequence': sequences,
        'match_chrm': ['chr1'] * len(sequences),
        'match_position': list(range(len(sequences))),
        'match_strand': ['+'] * len(sequences),
        'specificity': specificities,
    })


def test_apply_filters_removes_tttt():
    df = _make_guide_df(['ACGTTTTTACGTACGTACGT', 'ACGTACGTACGTACGTACGT'])
    result = apply_filters(df, specificity_thresh=0.2)
    assert len(result) == 1
    assert result.iloc[0]['sequence'] == 'ACGTACGTACGTACGTACGT'


def test_apply_filters_removes_bstxi():
    # BstXI site: CCACCTTGTTG (11nt) — embed at start of a 20nt spacer
    df = _make_guide_df(['CCACCTTGTTGACGTACGTA', 'ACGTACGTACGTACGTACGT'])
    result = apply_filters(df, specificity_thresh=0.2)
    assert len(result) == 1
    assert result.iloc[0]['sequence'] == 'ACGTACGTACGTACGTACGT'


def test_apply_filters_removes_bpi1102i():
    # Bpi1102I site: GTTTAAGAGCTAAGCTGG (18nt) — full site won't fit in 20nt,
    # but a partial match that is present should be caught.
    # Use a site that fits: first 18nt = GTTTAAGAGCTAAGCTGG + 2nt padding
    df = _make_guide_df(['GTTTAAGAGCTAAGCTGGAC', 'ACGTACGTACGTACGTACGT'])
    result = apply_filters(df, specificity_thresh=0.2)
    assert len(result) == 1
    assert result.iloc[0]['sequence'] == 'ACGTACGTACGTACGTACGT'


def test_apply_filters_removes_low_specificity():
    df = _make_guide_df(
        ['ACGTACGTACGTACGTACGT', 'TGCATGCATGCATGCATGCA'],
        specificities=[0.1, 0.9],
    )
    result = apply_filters(df, specificity_thresh=0.2)
    assert len(result) == 1
    assert result.iloc[0]['sequence'] == 'TGCATGCATGCATGCATGCA'


def test_apply_filters_removes_nan_specificity():
    df = _make_guide_df(['ACGTACGTACGTACGTACGT'], specificities=[float('nan')])
    result = apply_filters(df, specificity_thresh=0.2)
    assert len(result) == 0
```

- [ ] **Step 2: Run tests to confirm they fail**

```bash
uv run pytest testing/test_tiling.py -k "apply_filters" -v
```
Expected: FAIL (ImportError or NotImplementedError)

- [ ] **Step 3: Implement apply_filters in pyguide/tiling.py**

Add after `run_guidescan`:
```python
_BSTXI = "CCACCTTGTTG"
_BPI1102I = "GTTTAAGAGCTAAGCTGG"


def apply_filters(df: pd.DataFrame, specificity_thresh: float) -> pd.DataFrame:
    """Apply all sequence and specificity filters. Returns filtered DataFrame."""
    mask = (
        ~df['sequence'].str.contains('TTTT', na=False)
        & ~df['sequence'].str.contains(_BSTXI, na=False)
        & ~df['sequence'].str.contains(_BPI1102I, na=False)
        & df['specificity'].notna()
        & (df['specificity'] > specificity_thresh)
    )
    return df[mask].reset_index(drop=True)
```

- [ ] **Step 4: Run tests to confirm they pass**

```bash
uv run pytest testing/test_tiling.py -k "apply_filters" -v
```
Expected: all PASS

- [ ] **Step 5: Commit**

```bash
git add pyguide/tiling.py testing/test_tiling.py
git commit -m "feat: add apply_filters (TTTT, restriction sites, specificity)"
```

---

## Task 6: greedy_hamming_select

**Files:**
- Modify: `pyguide/tiling.py`
- Modify: `testing/test_tiling.py`

- [ ] **Step 1: Write the failing tests**

Add to `testing/test_tiling.py`:
```python
from pyguide.tiling import hamming_distance, greedy_hamming_select


def test_hamming_distance_identical():
    assert hamming_distance("ACGT", "ACGT") == 0


def test_hamming_distance_all_different():
    assert hamming_distance("AAAA", "TTTT") == 4


def test_hamming_distance_one_mismatch():
    assert hamming_distance("ACGT", "ACGG") == 1


def test_greedy_hamming_select_no_pair_within_threshold():
    # Build a set of guides where pairs are >4 apart
    seqs = [
        'AAAAAAAAAAAAAAAAAAAA',  # specificity 0.9
        'AAAAAAAAAAAAAAAACCCC',  # hamming 4 from above — should be EXCLUDED
        'TTTTTTTTTTTTTTTTTTTT',  # hamming 20 from first — INCLUDED
    ]
    df = _make_guide_df(seqs, specificities=[0.9, 0.8, 0.7])
    result = greedy_hamming_select(df, min_hamming=4, guides_per_region=None)
    seqs_out = list(result['sequence'])
    for i in range(len(seqs_out)):
        for j in range(i + 1, len(seqs_out)):
            assert hamming_distance(seqs_out[i], seqs_out[j]) > 4, \
                f"Pair {i},{j} has hamming <= 4"


def test_greedy_hamming_select_priority():
    # Highest-specificity guide must always be in the result
    seqs = [
        'AAAAAAAAAAAAAAAAAAAA',  # specificity 0.5 — lower
        'TTTTTTTTTTTTTTTTTTTT',  # specificity 0.9 — highest, hamming=20 from first
    ]
    df = _make_guide_df(seqs, specificities=[0.5, 0.9])
    result = greedy_hamming_select(df, min_hamming=4, guides_per_region=None)
    assert 'TTTTTTTTTTTTTTTTTTTT' in list(result['sequence'])


def test_greedy_hamming_select_cap_even_spread():
    # 6 guides — all pairwise hamming > 4 (pure homopolymers + alternating)
    seqs = [
        'AAAAAAAAAAAAAAAAAAAA',  # 20 A's
        'TTTTTTTTTTTTTTTTTTTT',  # 20 T's, hamming=20 from #0
        'GGGGGGGGGGGGGGGGGGGG',  # 20 G's, hamming=20 from both
        'CCCCCCCCCCCCCCCCCCCC',  # 20 C's, hamming=20 from all above
        'ACACACACACACACACACAC',  # alternating AC, hamming > 4 from all above
        'TGTGTGTGTGTGTGTGTGTG',  # alternating TG, hamming > 4 from all above
    ]
    positions = [0, 10, 20, 30, 40, 50]
    specificities = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4]
    df = pd.DataFrame({
        'sequence': seqs,
        'match_chrm': ['chr1'] * 6,
        'match_position': positions,
        'match_strand': ['+'] * 6,
        'specificity': specificities,
    })
    result = greedy_hamming_select(df, min_hamming=4, guides_per_region=3)
    assert len(result) == 3
    pos_out = sorted(result['match_position'].tolist())
    # Even spread: should include first (pos 0) and last (pos 50)
    assert pos_out[0] == 0
    assert pos_out[-1] == 50


def test_greedy_hamming_select_cap_larger_than_available():
    # If fewer guides available than cap, return all
    seqs = ['AAAAAAAAAAAAAAAAAAAA', 'TTTTTTTTTTTTTTTTTTTT']
    df = _make_guide_df(seqs, specificities=[0.9, 0.8])
    result = greedy_hamming_select(df, min_hamming=4, guides_per_region=10)
    assert len(result) == 2
```

- [ ] **Step 2: Run tests to confirm they fail**

```bash
uv run pytest testing/test_tiling.py -k "hamming" -v
```
Expected: FAIL (ImportError)

- [ ] **Step 3: Implement hamming_distance and greedy_hamming_select in pyguide/tiling.py**

Add after `apply_filters`:
```python
def hamming_distance(s1: str, s2: str) -> int:
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def greedy_hamming_select(
    df: pd.DataFrame,
    min_hamming: int,
    guides_per_region: Optional[int],
) -> pd.DataFrame:
    """
    Greedily select guides sorted by descending specificity such that
    no two selected guides have pairwise Hamming distance <= min_hamming.

    If guides_per_region is set, subsample the selected guides at evenly
    spaced genomic positions (maximizes tiling coverage).
    """
    if df.empty:
        return df

    df_sorted = df.sort_values('specificity', ascending=False).reset_index(drop=True)
    selected_rows = []
    selected_seqs = []

    for _, row in df_sorted.iterrows():
        seq = row['sequence']
        if all(hamming_distance(seq, s) > min_hamming for s in selected_seqs):
            selected_rows.append(row)
            selected_seqs.append(seq)

    if not selected_rows:
        return df.iloc[0:0]  # empty with same columns

    result = pd.DataFrame(selected_rows).reset_index(drop=True)

    if guides_per_region is not None and len(result) > guides_per_region:
        result_by_pos = result.sort_values('match_position').reset_index(drop=True)
        n = len(result_by_pos)
        N = guides_per_region
        indices = [round(i * (n - 1) / (N - 1)) for i in range(N)]
        result = result_by_pos.iloc[indices].reset_index(drop=True)

    return result
```

- [ ] **Step 4: Run tests to confirm they pass**

```bash
uv run pytest testing/test_tiling.py -k "hamming" -v
```
Expected: all PASS

- [ ] **Step 5: Commit**

```bash
git add pyguide/tiling.py testing/test_tiling.py
git commit -m "feat: add hamming_distance and greedy_hamming_select with positional cap"
```

---

## Task 7: write_sequence_file

**Files:**
- Modify: `pyguide/tiling.py`
- Modify: `testing/test_tiling.py`

- [ ] **Step 1: Write the failing tests**

Add to `testing/test_tiling.py`:
```python
from pyguide.tiling import write_sequence_file


def test_write_sequence_file_format(tmp_path):
    df = pd.DataFrame({
        'sequence': ['ACGTACGTACGTACGTACGT', 'TGCATGCATGCATGCATGCA'],
        'match_chrm': ['chr1', 'chr1'],
        'match_position': [12350, 12370],
        'match_strand': ['+', '-'],
        'specificity': [0.9, 0.8],
    })
    out = str(tmp_path / "out.txt")
    write_sequence_file(df, out)
    lines = open(out).read().strip().split('\n')
    assert len(lines) == 2
    parts0 = lines[0].split('\t')
    assert len(parts0) == 2
    assert parts0[0] == 'chr1:12350_fwd'
    assert parts0[1] == 'ACGTACGTACGTACGTACGT'
    parts1 = lines[1].split('\t')
    assert parts1[0] == 'chr1:12370_rev'
    assert parts1[1] == 'TGCATGCATGCATGCATGCA'


def test_write_sequence_file_empty(tmp_path):
    df = pd.DataFrame(columns=['sequence', 'match_chrm', 'match_position', 'match_strand'])
    out = str(tmp_path / "out.txt")
    write_sequence_file(df, out)
    assert open(out).read() == ""
```

- [ ] **Step 2: Run tests to confirm they fail**

```bash
uv run pytest testing/test_tiling.py -k "write_sequence_file" -v
```
Expected: FAIL (ImportError)

- [ ] **Step 3: Implement write_sequence_file in pyguide/tiling.py**

Add after `greedy_hamming_select`:
```python
def write_sequence_file(df: pd.DataFrame, output_path: str) -> None:
    """Write tab-delimited name<TAB>sequence file. Names encode chrom:position_strand."""
    with open(output_path, 'w') as f:
        for _, row in df.iterrows():
            strand_label = 'fwd' if row['match_strand'] == '+' else 'rev'
            name = f"{row['match_chrm']}:{row['match_position']}_{strand_label}"
            f.write(f"{name}\t{row['sequence']}\n")
```

- [ ] **Step 4: Run tests to confirm they pass**

```bash
uv run pytest testing/test_tiling.py -k "write_sequence_file" -v
```
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add pyguide/tiling.py testing/test_tiling.py
git commit -m "feat: add write_sequence_file with chrom:pos_strand naming"
```

---

## Task 8: write_bigwig

**Files:**
- Modify: `pyguide/tiling.py`
- Modify: `testing/test_tiling.py`

- [ ] **Step 1: Write the failing tests**

Add to `testing/test_tiling.py`:
```python
from pyguide.tiling import write_bigwig


def test_write_bigwig_missing_pybigwig(tmp_path):
    df = pd.DataFrame({
        'sequence': ['ACGTACGTACGTACGTACGT'],
        'match_chrm': ['chr1'],
        'match_position': [12350],
        'match_strand': ['+'],
        'specificity': [0.9],
    })
    chrom_sizes = str(tmp_path / "chrom.sizes")
    with open(chrom_sizes, 'w') as f:
        f.write("chr1\t248956422\n")
    out = str(tmp_path / "out.bw")

    with patch.dict('sys.modules', {'pyBigWig': None}):
        with pytest.raises(ImportError, match="uv sync --extra tiling"):
            write_bigwig(df, out, chrom_sizes)


def test_write_bigwig_creates_file(tmp_path):
    """Only runs if pyBigWig is installed; skips otherwise."""
    pytest.importorskip("pyBigWig")
    import pyBigWig

    df = pd.DataFrame({
        'sequence': ['ACGTACGTACGTACGTACGT'],
        'match_chrm': ['chr1'],
        'match_position': [100],  # 1-based
        'match_strand': ['+'],
        'specificity': [0.9],
    })
    chrom_sizes = str(tmp_path / "chrom.sizes")
    with open(chrom_sizes, 'w') as f:
        f.write("chr1\t248956422\n")
    out = str(tmp_path / "out.bw")
    write_bigwig(df, out, chrom_sizes)
    assert os.path.exists(out)

    bw = pyBigWig.open(out)
    # Guide at 1-based position 100 covers 0-based positions 99..118
    vals = bw.values("chr1", 99, 119)
    bw.close()
    assert all(v == 1.0 for v in vals)
```

- [ ] **Step 2: Run tests to confirm they fail**

```bash
uv run pytest testing/test_tiling.py -k "bigwig" -v
```
Expected: FAIL (ImportError)

- [ ] **Step 3: Implement write_bigwig in pyguide/tiling.py**

Add after `write_sequence_file`:
```python
def write_bigwig(df: pd.DataFrame, output_path: str, chrom_sizes_path: str) -> None:
    """
    Write a BigWig coverage file. Each guide's 20nt spacer contributes 1 per base.
    Requires pyBigWig: uv sync --extra tiling
    """
    try:
        import pyBigWig
    except ImportError:
        raise ImportError(
            "pyBigWig is required for BigWig output. "
            "Install it with: uv sync --extra tiling"
        )

    # Load chrom sizes
    chrom_sizes = {}
    with open(chrom_sizes_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                chrom_sizes[parts[0]] = int(parts[1])

    # Accumulate coverage per chromosome
    from collections import defaultdict
    coverage: dict[str, dict[int, int]] = defaultdict(lambda: defaultdict(int))

    for _, row in df.iterrows():
        chrom = row['match_chrm']
        pos = int(row['match_position'])  # 1-based
        start_0 = pos - 1  # convert to 0-based
        chrom_len = chrom_sizes.get(chrom, 0)
        for base in range(start_0, min(start_0 + 20, chrom_len)):
            coverage[chrom][base] += 1

    bw = pyBigWig.open(output_path, 'w')
    # BigWig header must list all chroms that will appear
    present_chroms = [(c, chrom_sizes[c]) for c in chrom_sizes if c in coverage]
    bw.addHeader(present_chroms)

    for chrom, pos_dict in sorted(coverage.items()):
        positions = sorted(pos_dict.keys())
        # Write as a series of individual intervals
        starts = positions
        ends = [p + 1 for p in positions]
        values = [float(pos_dict[p]) for p in positions]
        bw.addEntries(
            [chrom] * len(starts),
            starts,
            ends=ends,
            values=values,
        )

    bw.close()
    print(f"BigWig written: {output_path}")
```

- [ ] **Step 4: Run tests to confirm they pass**

```bash
uv run pytest testing/test_tiling.py -k "bigwig" -v
```
Expected: PASS (test_write_bigwig_missing_pybigwig always runs; test_write_bigwig_creates_file skips if pyBigWig not installed)

- [ ] **Step 5: Commit**

```bash
git add pyguide/tiling.py testing/test_tiling.py
git commit -m "feat: add write_bigwig with pyBigWig optional dependency"
```

---

## Task 9: main() CLI entry point

**Files:**
- Modify: `pyguide/tiling.py`
- Modify: `testing/test_tiling.py`

- [ ] **Step 1: Write the failing integration test**

Add to `testing/test_tiling.py`:
```python
from pyguide.tiling import main as tiling_main


def test_main_requires_coordinates_or_file(tmp_path, capsys):
    with pytest.raises(SystemExit):
        tiling_main(["--index", "/fake/index", "--output", str(tmp_path / "out.txt")])


def test_main_coordinates_and_file_mutually_exclusive(tmp_path, capsys):
    coords_file = tmp_path / "coords.txt"
    coords_file.write_text("chr1:12345-12395\n")
    with pytest.raises(SystemExit):
        tiling_main([
            "--index", "/fake/index",
            "--output", str(tmp_path / "out.txt"),
            "--coordinates", "chr1:12345-12395",
            "--coordinates_file", str(coords_file),
        ])
```

- [ ] **Step 2: Run tests to confirm they fail**

```bash
uv run pytest testing/test_tiling.py -k "test_main" -v
```
Expected: FAIL

- [ ] **Step 3: Implement main() in pyguide/tiling.py**

Replace the existing stub `main()`:
```python
def main(raw_args=None) -> None:
    parser = ArgumentParser(
        description="Generate a tiling gRNA library from chromosomal coordinates using GuideScan2."
    )

    coord_group = parser.add_mutually_exclusive_group(required=True)
    coord_group.add_argument(
        "--coordinates",
        type=str,
        help="Single coordinate region, e.g. chr1:12345-12395 (1-based, closed).",
    )
    coord_group.add_argument(
        "--coordinates_file",
        type=str,
        help="File with one coordinate region per line.",
    )

    parser.add_argument("--index", type=str, required=True,
                        help="Path to GuideScan2 hg38 index file.")
    parser.add_argument("--output", type=str, required=True,
                        help="Output tab-delimited sequence file (name<TAB>sequence).")
    parser.add_argument("--specificity", type=float, default=0.2,
                        help="Minimum GuideScan2 specificity score (default: 0.2).")
    parser.add_argument("--hamming", type=int, default=4,
                        help="Minimum pairwise Hamming distance between guides (default: 4).")
    parser.add_argument("--guides_per_region", type=int, default=None,
                        help="Max guides per region, evenly spread (optional).")
    parser.add_argument("--bigwig", action="store_true",
                        help="Also write a BigWig coverage file alongside --output.")

    args = parser.parse_args(raw_args)

    if not os.path.exists(args.index):
        print(f"Error: index file not found: {args.index}", file=sys.stderr)
        sys.exit(1)

    # Collect coordinate strings
    if args.coordinates:
        coord_strings = [args.coordinates]
    else:
        with open(args.coordinates_file) as f:
            coord_strings = [line.strip() for line in f if line.strip()]

    all_guides: list[pd.DataFrame] = []

    for coord_str in coord_strings:
        try:
            chrom, start, end = parse_coordinates(coord_str)
        except ValueError as e:
            print(f"Warning: skipping invalid coordinate '{coord_str}': {e}", file=sys.stderr)
            continue

        print(f"Processing {coord_str}...")

        try:
            sequence = fetch_sequence_ucsc(chrom, start, end)
        except RuntimeError as e:
            print(f"Warning: skipping {coord_str}: {e}", file=sys.stderr)
            continue

        raw_guides = find_ngg_guides(chrom, start, sequence)
        if not raw_guides:
            print(f"Warning: no NGG guides found in {coord_str}", file=sys.stderr)
            continue

        try:
            guides_df = run_guidescan(raw_guides, args.index)
        except RuntimeError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)

        filtered = apply_filters(guides_df, specificity_thresh=args.specificity)
        if filtered.empty:
            print(f"Warning: no guides passed filters for {coord_str}", file=sys.stderr)
            continue

        selected = greedy_hamming_select(filtered, args.hamming, args.guides_per_region)
        if selected.empty:
            print(f"Warning: no guides survived Hamming selection for {coord_str}", file=sys.stderr)
            continue

        print(f"  {len(selected)} guides selected for {coord_str}")
        all_guides.append(selected)

    if not all_guides:
        print("No guides selected for any region. Output file not written.", file=sys.stderr)
        sys.exit(1)

    combined = pd.concat(all_guides, ignore_index=True)
    write_sequence_file(combined, args.output)
    print(f"Sequence file written: {args.output}")

    if args.bigwig:
        bw_path = os.path.splitext(args.output)[0] + ".bw"
        here = os.path.dirname(os.path.realpath(__file__))
        chrom_sizes_path = os.path.join(here, 'data', 'hg38.chrom.sizes')
        write_bigwig(combined, bw_path, chrom_sizes_path)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests to confirm they pass**

```bash
uv run pytest testing/test_tiling.py -k "test_main" -v
```
Expected: PASS

- [ ] **Step 5: Run the full test suite to check for regressions**

```bash
uv run pytest testing/test_tiling.py -v
```
Expected: all tests PASS

- [ ] **Step 6: Commit**

```bash
git add pyguide/tiling.py testing/test_tiling.py
git commit -m "feat: add main() CLI entry point for pyguide-tiling"
```

---

## Task 10: Download script

**Files:**
- Create: `scripts/download_guidescan_data.sh`

- [ ] **Step 1: Create the scripts directory and download script**

```bash
mkdir -p scripts
```

Write `scripts/download_guidescan_data.sh`:
```bash
#!/usr/bin/env bash
# Download GuideScan2 hg38 index (~2.2 GB) and hg38 chrom sizes.
# Usage: bash scripts/download_guidescan_data.sh [--data-dir /path/to/dir]

set -euo pipefail

DATA_DIR="${HOME}/.pyguide_data"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1" >&2
            echo "Usage: $0 [--data-dir /path/to/dir]" >&2
            exit 1
            ;;
    esac
done

mkdir -p "${DATA_DIR}"

INDEX_URL="https://guidescan.com/indices/hg38.zip"
INDEX_ZIP="${DATA_DIR}/hg38.zip"
INDEX_FILE="${DATA_DIR}/hg38.index"

CHROM_SIZES_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"
CHROM_SIZES_DEST="$(dirname "$0")/../pyguide/data/hg38.chrom.sizes"

echo "=== GuideScan2 hg38 Data Download ==="
echo "Data directory: ${DATA_DIR}"
echo ""

# Download hg38 GuideScan2 index
if [[ -f "${INDEX_FILE}" ]]; then
    echo "[SKIP] hg38 index already exists: ${INDEX_FILE}"
else
    echo "[DOWNLOAD] hg38 GuideScan2 index (~2.2 GB) ..."
    echo "  Note: Using --insecure due to expired SSL cert on guidescan.com"
    curl --insecure -L --progress-bar -o "${INDEX_ZIP}" "${INDEX_URL}"
    echo "[EXTRACT] Extracting hg38.zip ..."
    unzip -o "${INDEX_ZIP}" -d "${DATA_DIR}"
    # Find the extracted index file (may have a different name inside zip)
    EXTRACTED=$(find "${DATA_DIR}" -name "*.index" -o -name "hg38*" ! -name "*.zip" 2>/dev/null | head -1)
    if [[ -z "${EXTRACTED}" ]]; then
        echo "ERROR: Could not find index file after extraction." >&2
        echo "Contents of ${DATA_DIR}:" >&2
        ls "${DATA_DIR}" >&2
        exit 1
    fi
    echo "[OK] Index extracted: ${EXTRACTED}"
fi

# Download hg38 chrom sizes (for BigWig output)
echo ""
DEST_DIR="$(dirname "${CHROM_SIZES_DEST}")"
mkdir -p "${DEST_DIR}"
if [[ -f "${CHROM_SIZES_DEST}" ]]; then
    echo "[SKIP] hg38.chrom.sizes already exists: ${CHROM_SIZES_DEST}"
else
    echo "[DOWNLOAD] hg38 chromosome sizes ..."
    curl -L --progress-bar -o "${CHROM_SIZES_DEST}" "${CHROM_SIZES_URL}"
    echo "[OK] chrom sizes written: ${CHROM_SIZES_DEST}"
fi

echo ""
echo "=== Done ==="
echo "Index location: ${DATA_DIR}"
echo "Pass to pyguide-tiling with:  --index ${DATA_DIR}/<index-filename>"
echo "Run 'ls ${DATA_DIR}' to see the extracted index filename."
```

- [ ] **Step 2: Make script executable**

```bash
chmod +x scripts/download_guidescan_data.sh
```

- [ ] **Step 3: Verify the script is syntactically valid**

```bash
bash -n scripts/download_guidescan_data.sh && echo "Syntax OK"
```
Expected: `Syntax OK`

- [ ] **Step 4: Commit**

```bash
git add scripts/download_guidescan_data.sh
git commit -m "feat: add download script for GuideScan2 hg38 index"
```

---

## Task 11: README update

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Add GuideScan2 tiling section to README.md**

Open `README.md` and add a new section after the existing `pyguide-check-seq` usage section. Add the following:

````markdown
## Genomic Tiling Library Usage

### Prerequisites

**GuideScan2** must be installed separately (not available on PyPI):
```bash
conda install -c bioconda guidescan2
```

**hg38 GuideScan2 index** (~2.2 GB) — download with the provided script:
```bash
bash scripts/download_guidescan_data.sh
# Default saves to ~/.pyguide_data/
# Use --data-dir to choose a custom location
```

**BigWig output** (optional) — install pyBigWig:
```bash
uv sync --extra tiling
```

---

### Flags for pyguide-tiling

**Required (one of):**
- **--coordinates**: Single coordinate region, e.g. `chr1:12345-12395` (1-based, closed)
- **--coordinates_file**: Path to a file with one coordinate region per line

**Required:**
- **--index**: Path to the GuideScan2 hg38 index file downloaded above
- **--output**: Path for the output tab-delimited sequence file

**Optional:**
- **--specificity**: Minimum GuideScan2 specificity score (default: 0.2). Higher = more specific guides only.
- **--hamming**: Minimum pairwise Hamming distance between any two selected guides (default: 4). Prevents near-identical oligos in the pool.
- **--guides_per_region**: Maximum number of guides per region, evenly spread across positions to maximize tiling coverage (optional).
- **--bigwig**: Also write a BigWig coverage track (`.bw` file) alongside the output, viewable in IGV or UCSC Genome Browser. Requires `uv sync --extra tiling`.

---

### Example use cases

**(1) Tile a single promoter region:**
```bash
pyguide-tiling \
  --index ~/.pyguide_data/hg38.index \
  --coordinates chr1:12345-12595 \
  --output tiling_guides.txt
```

**(2) Tile multiple regions from a file:**
```bash
# regions.txt contains one region per line:
# chr1:12345-12595
# chr7:117548601-117548801

pyguide-tiling \
  --index ~/.pyguide_data/hg38.index \
  --coordinates_file regions.txt \
  --output tiling_guides.txt \
  --guides_per_region 20
```

**(3) Tile with BigWig coverage output:**
```bash
pyguide-tiling \
  --index ~/.pyguide_data/hg38.index \
  --coordinates chr7:117548601-117548801 \
  --output tiling_guides.txt \
  --bigwig
# Produces: tiling_guides.txt and tiling_guides.bw
```

### pyguide-tiling output

A tab-delimited file with two columns: `name<TAB>20nt_spacer`. Guide names encode genomic position and strand (e.g. `chr1:12350_fwd`, `chr1:12370_rev`).

This file feeds directly into the pooled-seq ordering pipeline:
```bash
# Step 1: generate tiling guides
pyguide-tiling --index ~/.pyguide_data/hg38.index \
  --coordinates chr1:12345-12595 --output tiling_guides.txt

# Step 2: collate (assigns library primers)
pyguide-collate-seq --sequence_files tiling_guides.txt

# Step 3: order the pooled oligo library
pyguide-order --wishlist_file collated_seq_wishlist_YY_MM_DD.txt \
  --name Name --ai i --guides_per_gene 1 --order_format pooled-seq
```
````

- [ ] **Step 2: Also add `tiling` to the Installation section**

In the existing Installation section of `README.md`, add after the `**Everything**` block:

```markdown
**Tiling library + BigWig output** (includes pyBigWig for genome browser tracks):
```bash
uv sync --extra tiling
```
```

- [ ] **Step 3: Commit**

```bash
git add README.md
git commit -m "docs: add pyguide-tiling usage, GuideScan2 setup, and tiling installation to README"
```

---

## Self-Review

### Spec Coverage Check

| Spec requirement | Task |
|---|---|
| `pyguide-tiling` CLI command | Task 1, 9 |
| `--coordinates` and `--coordinates_file` (mutually exclusive) | Task 9 |
| GuideScan2 subprocess via `guidescan enumerate` | Task 4 |
| Runtime check for `guidescan` on PATH | Task 4 |
| UCSC REST API for sequence fetching | Task 3 |
| NGG PAM scanning (both strands) | Task 3 |
| TTTT motif filter | Task 5 |
| BstXI restriction site filter | Task 5 |
| Bpi1102I restriction site filter | Task 5 |
| Specificity > threshold filter | Task 5 |
| Greedy pairwise Hamming > 4 selection | Task 6 |
| Selection sorted by descending specificity | Task 6 |
| Uniform positional subsample cap | Task 6 |
| `name<TAB>sequence` output with chrom:pos_strand naming | Task 7 |
| BigWig coverage output (optional, `--bigwig`) | Task 8 |
| pyBigWig optional dep, clear error if missing | Task 8 |
| `tiling = ["pyBigWig"]` optional dep in pyproject.toml | Task 1 |
| Bundle decode_database.py + CFD pickles (MIT) | Task 2 |
| Bundle hg38.chrom.sizes | Task 2 |
| Download script for hg38 index + chrom sizes | Task 10 |
| README: GuideScan2 setup + tiling usage | Task 11 |
| Per-region warning (not abort) when no guides survive | Task 9 |
| RuntimeError with clear message for missing guidescan | Task 4 |

All spec requirements are covered.
