import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.request
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


def write_sequence_file(df: pd.DataFrame, output_path: str) -> None:
    """Write tab-delimited name<TAB>sequence file. Names encode chrom:position_strand."""
    with open(output_path, 'w') as f:
        for _, row in df.iterrows():
            strand_label = 'fwd' if row['match_strand'] == '+' else 'rev'
            name = f"{row['match_chrm']}:{row['match_position']}_{strand_label}"
            f.write(f"{name}\t{row['sequence']}\n")


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


def main(raw_args=None):
    pass
