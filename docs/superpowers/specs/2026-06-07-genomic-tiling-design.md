# Genomic Tiling Library Feature — Design Spec

**Date:** 2026-06-07
**Branch:** `genomic_tiling`

---

## Overview

Add a `pyguide-tiling` CLI command that accepts chromosomal coordinates, runs GuideScan2 to enumerate all SpCas9 NGG-PAM guides in the region, applies quality and specificity filters, performs greedy pairwise hamming distance selection, and outputs a tab-delimited sequence file compatible with the existing `pyguide-collate-seq` → `pyguide-order --order_format pooled-seq` pipeline.

Optionally generates a BigWig coverage track for visualization in genome browsers (IGV, UCSC).

---

## External Dependencies

### GuideScan2 (required, not on PyPI)

GuideScan2 must be installed separately and available on `PATH` as `guidescan`. It is a standalone bioinformatics CLI tool. See README for installation instructions (Conda recommended).

At runtime, `pyguide-tiling` checks for `guidescan` on PATH and fails fast with a clear message if not found.

### hg38 Genome Index (required, ~2.2 GB)

Downloaded via `scripts/download_guidescan_data.sh`. Stored in user-specified data directory (default: `~/.pyguide_data/`).

URL: `https://guidescan.com/indices/hg38.zip` (SSL cert may be expired — script uses `--insecure` flag)

### pyBigWig (optional Python dependency)

Required only for `--bigwig` output. Installed via `uv sync --extra tiling`. If `--bigwig` is passed without pyBigWig installed, the command exits with a clear install instruction.

### hg38 Chromosome Sizes (bundled)

`pyguide/data/hg38.chrom.sizes` — small text file (~3 KB), bundled in the repo. Also downloaded/refreshed by the download script.

### GuideScan2 CFD Scripts (bundled, MIT license)

Copied from https://github.com/pritykinlab/guidescan-cli/tree/master/scripts into `pyguide/tiling_scripts/`:
- `decode_database.py`
- `cfd/mismatch_score.pkl`
- `cfd/pam_scores.pkl`

Bundled for reference and future extension. The main pipeline uses GuideScan2's `specificity` CSV column directly (which is the CFD-based aggregate score: `1 / (1 + CFD_sum)`).

---

## Architecture & Data Flow

```
pyguide-tiling
  --index /path/to/hg38.index
  --coordinates chr1:12345-12395        (or --coordinates_file)
  --output /path/to/sequences.txt
  [--specificity 0.2]
  [--hamming 4]
  [--guides_per_region N]
  [--bigwig]
        │
        ▼ (one subprocess call per region)
guidescan enumerate --index ... --pam NGG
        │  CSV output: id, sequence, chromosome, position, sense, specificity, ...
        ▼
Filter pipeline (pandas, in-memory):
  1. Remove guides with TTTT motif in spacer
  2. Remove guides with BstXI site ("CCACCTTGTTG") in spacer
  3. Remove guides with Bpi1102I site ("GTTTAAGAGCTAAGCTGG") in spacer
  4. Keep guides with specificity > threshold (default 0.2)
  5. Greedy pairwise hamming > 4 selection (sorted by descending specificity)
  6. If --guides_per_region N: uniform positional subsample → N evenly spread guides
        │
        ▼
Tab-delimited output: name<TAB>20nt_spacer
  Name format: chr1:12350_fwd  or  chr1:12350_rev
        │
        ▼  (optional)
BigWig coverage file (per-base coverage, each 20nt spacer contributes 1 per base)
        │
        ▼
pyguide-collate-seq → pyguide-order --order_format pooled-seq
```

**No NGG PAM check in Python** — GuideScan2 enforces NGG at enumeration time via `--pam NGG`.

---

## Module Design: `pyguide/tiling.py`

### Functions

**`parse_coordinates(coord_str: str) → tuple[str, int, int]`**
- Input: `"chr1:12345-12395"`
- Output: `("chr1", 12345, 12395)`
- Raises `ValueError` with clear message on malformed input

**`run_guidescan(chrom: str, start: int, end: int, index_path: str) → pd.DataFrame`**
- Checks `guidescan` is on PATH; raises `RuntimeError` with install message if not
- Shells out: `guidescan enumerate --index <index> --pam NGG --chromosome <chrom> --start <start> --end <end>`
- Returns parsed CSV as DataFrame
- Raises `RuntimeError` if subprocess fails

**`apply_filters(df: pd.DataFrame, specificity_thresh: float) → pd.DataFrame`**
- Removes rows where spacer contains `TTTT`
- Removes rows where spacer contains `CCACCTTGTTG` (BstXI)
- Removes rows where spacer contains `GTTTAAGAGCTAAGCTGG` (Bpi1102I)
- Removes rows where `specificity <= specificity_thresh`
- Returns filtered DataFrame

**`greedy_hamming_select(df: pd.DataFrame, min_hamming: int, guides_per_region: int | None) → pd.DataFrame`**
- Sorts by `specificity` descending
- Greedy loop: add guide only if hamming distance > `min_hamming` from every already-selected guide
- If `guides_per_region` is set: sort selected guides by genomic position, uniformly subsample N (every k-th, k = len(selected) // N) to maximize positional spread
- Returns selected DataFrame

**`write_sequence_file(df: pd.DataFrame, output_path: str) → None`**
- Writes tab-delimited `name\tsequence` per line
- Name format: `{chrom}:{position}_fwd` or `{chrom}:{position}_rev`

**`write_bigwig(df: pd.DataFrame, output_path: str, chrom_sizes_path: str) → None`**
- Imports `pyBigWig`; raises `ImportError` with install message if missing
- Each guide covers its 20nt spacer positions (position to position+19)
- Writes per-base integer coverage track

**`main() → None`**
- argparse entry point
- Reads `--coordinates` or `--coordinates_file` (mutually exclusive, one required)
- Loops over regions, calls run_guidescan → apply_filters → greedy_hamming_select
- Concatenates results across regions, writes sequence file
- If `--bigwig`: calls write_bigwig with bundled `pyguide/data/hg38.chrom.sizes`
- Per-region warning (not abort) if no guides survive filtering

---

## CLI Interface

```
pyguide-tiling
  --index PATH              Path to GuideScan2 hg38 index file (required)
  --coordinates STR         Single coordinate region, e.g. chr1:12345-12395
  --coordinates_file PATH   File with one coordinate region per line
  --output PATH             Output tab-delimited sequence file (required)
  --specificity FLOAT       Minimum specificity score, default 0.2
  --hamming INT             Minimum pairwise hamming distance, default 4
  --guides_per_region INT   Max guides per region, evenly spread (optional)
  --bigwig                  Also write BigWig coverage file alongside --output
```

`--coordinates` and `--coordinates_file` are mutually exclusive; exactly one is required.

---

## New Files

| Path | Purpose |
|---|---|
| `pyguide/tiling.py` | Main module + `main()` entry point |
| `pyguide/tiling_scripts/decode_database.py` | Bundled from guidescan-cli (MIT) |
| `pyguide/tiling_scripts/cfd/mismatch_score.pkl` | Bundled from guidescan-cli (MIT) |
| `pyguide/tiling_scripts/cfd/pam_scores.pkl` | Bundled from guidescan-cli (MIT) |
| `pyguide/data/hg38.chrom.sizes` | Bundled hg38 chromosome sizes for BigWig |
| `scripts/download_guidescan_data.sh` | Downloads hg38 index + chrom sizes |
| `testing/test_tiling.py` | Unit tests |

## Modified Files

| Path | Change |
|---|---|
| `pyproject.toml` | Add `pyguide-tiling` entry point; add `tiling = ["pyBigWig"]` optional dep |
| `README.md` | Add GuideScan2 setup section, tiling usage, BigWig output docs |

---

## Filter Pipeline Detail

### 1. TTTT Motif
Guides containing four or more consecutive T's cause early transcriptional termination from the U6 Pol III promoter driving the guide cassette. Any spacer matching `.*TTTT.*` is removed.

### 2. Restriction Site Exclusion
Spacers must not contain the cloning restriction sites added by the pMK1334 protocol:
- BstXI: `CCACCTTGTTG`
- Bpi1102I: `GTTTAAGAGCTAAGCTGG`

### 3. Specificity Filter
GuideScan2 reports `specificity = 1 / (1 + CFD_sum)` where CFD_sum is the sum of CFD scores across all off-targets. Higher = more specific. Default threshold: `specificity > 0.2`.

### 4. Greedy Pairwise Hamming Selection
To avoid near-identical oligos in the pool that could cross-hybridize or confuse sequencing deconvolution:
- Sort candidates by `specificity` descending (best guide first)
- Greedily add a guide only if its Hamming distance to every already-selected guide exceeds `min_hamming` (default 4)
- Result: maximally specific, non-redundant set

### 5. Uniform Positional Subsampling (optional cap)
If `--guides_per_region N` is set, after greedy selection:
- Sort selected guides by genomic position
- Pick N guides at evenly spaced indices: `indices = [round(i * (len(selected)-1) / (N-1)) for i in range(N)]`
- Ensures even coverage across the region rather than clustering at the highest-specificity end
- If `len(selected) <= N`, returns all selected guides

---

## Error Handling

| Condition | Behavior |
|---|---|
| `guidescan` not on PATH | `RuntimeError`: clear message + README link |
| Index file not found | `FileNotFoundError`: path shown |
| Malformed coordinate string | `ValueError`: expected format shown |
| No guides after filtering (per region) | Warning printed, region skipped, execution continues |
| `--bigwig` without pyBigWig | `ImportError`: `uv sync --extra tiling` instruction |
| guidescan subprocess fails | `RuntimeError`: stderr captured and shown |

---

## Testing

**`testing/test_tiling.py`**

- `test_parse_coordinates_valid` — standard input parses correctly
- `test_parse_coordinates_invalid` — malformed strings raise ValueError
- `test_apply_filters_tttt` — guides with TTTT removed
- `test_apply_filters_bstxi` — guides with BstXI site removed
- `test_apply_filters_bpi1102i` — guides with Bpi1102I site removed
- `test_apply_filters_specificity` — guides below threshold removed
- `test_greedy_hamming_select_pairwise` — no pair in output has hamming ≤ 4
- `test_greedy_hamming_select_priority` — highest-specificity guide always selected first
- `test_greedy_hamming_select_cap` — N guides are evenly spread across region positions
- `test_write_sequence_file_format` — tab-delimited, correct name encoding (chrom:pos_strand)
- `test_run_guidescan_missing` — mock PATH empty, verify RuntimeError with install message

---

## Download Script: `scripts/download_guidescan_data.sh`

```
Usage: bash scripts/download_guidescan_data.sh [--data-dir /path/to/dir]

Downloads:
  - hg38 GuideScan2 index (~2.2 GB) from guidescan.com/indices/hg38.zip
  - hg38.chrom.sizes and copies to pyguide/data/hg38.chrom.sizes

Behavior:
  - Default data directory: ~/.pyguide_data/
  - Skips download if file already exists (idempotent)
  - Uses curl --insecure to bypass expired SSL cert on guidescan.com
  - Extracts zip, confirms index file present
```

---

## README Additions

- **Installation**: new `tiling` extra — `uv sync --extra tiling`
- **GuideScan2 setup**: install instructions (Conda), run `scripts/download_guidescan_data.sh`
- **`pyguide-tiling` usage**: example commands for single coordinate, coordinate file, with `--bigwig`
- **Pipeline**: show full tiling → collate-seq → order-pooled-seq workflow
