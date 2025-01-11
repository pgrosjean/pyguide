import numpy as np
import pandas as pd
import os
from datetime import datetime
from typing import List, Tuple
from sys import platform
from argparse import ArgumentParser
from pyguide.guide import get_unique_filename, read_wishlist
from pyguide.pool import get_primer_list


def generate_batch_retest_list(wishlist_file: str,
                               primer_file: str):
    """
    Generates a batch retesting input file based on guide IDs from a wishlist.

    Parameters
    ----------
    wishlist_file : str
        Path to the wishlist file containing guide IDs.
    primer_file : str
        Path to the primer file with primer pairs. If None, primers are assigned automatically.

        
    Returns
    -------
    None
    """
    # Read guide IDs from wishlist
    guide_ids, _ = read_wishlist(wishlist_file)  # Assumes only guide IDs are present
    assert len(guide_ids) > 0, "The wishlist file must contain guide IDs."

    # Assign primers
    if primer_file is None:
        # Use default primer file path
        file_path = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
        primer_file = os.path.join(file_path, '..', 'data', 'pooled_primers.txt')
    primer_list = get_primer_list(primer_file, 1)  # Only need one set of primers for batch retesting
    primer_tuple = primer_list[0]

    # Define output directory and file name
    base_dir = os.path.dirname(wishlist_file)
    now = datetime.now()
    date = now.strftime("%y_%m_%d")
    output_file_name = f"batch_retest_wishlist_{date}.txt"
    output_file_name = get_unique_filename(base_dir, output_file_name)
    output_file_path = os.path.join(base_dir, output_file_name)

    # Write batch retesting file
    primer_list = [primer_tuple] * len(guide_ids)
    with open(output_file_path, 'w') as outfile:
        lib_num = 0
        for guide_id, (left_primer, right_primer) in zip(guide_ids, primer_list):
            outfile.write(f"{guide_id}\t{left_primer}\t{right_primer}\t{lib_num}\n")
    
    print(f"Batch retest wishlist generated: {output_file_path}")


def main(raw_args=None):
    # Command-line argument parsing
    parser = ArgumentParser()
    parser.add_argument("--wishlist_file",
                        type=str,
                        required=True,
                        help="Path to the wishlist file containing guide IDs.")
    parser.add_argument("--primer_file",
                        type=str,
                        default=None,
                        help="Path to the primer file containing primer pairs. If not provided, default primers will be used.")

    # Parse arguments
    args = parser.parse_args(raw_args)

    # Run the batch retesting list generation
    generate_batch_retest_list(wishlist_file=args.wishlist_file,
                               primer_file=args.primer_file)


if __name__ == "__main__":
    main()
