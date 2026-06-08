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


def main(raw_args=None):
    pass
