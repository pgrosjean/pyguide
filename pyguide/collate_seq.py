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
