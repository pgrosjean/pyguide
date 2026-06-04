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
