import numpy as np
import pandas as pd
import os
from datetime import datetime
from typing import List, Tuple
from sys import platform
from argparse import ArgumentParser
from pyguide.guide import get_unique_filename, read_gene_list


####################
# Defining Functions
####################

## TODO: Expand the number of pooled primers that can be used

def generate_pooled_list(file_list: List[str],
                         primer_list: List[Tuple[str, str]]):
    """
    This function generates the input file for ordering a pooled library containing multiple sub-libraries.

    Parameters
    ----------
    file_list : List[str]
        List of gene of interest files.
    primer_list : List[Tuple[str, str]]
        List of primers for use in ordering pooled guides.

    Returns
    -------
    None

    """
    assert len(file_list) == len(primer_list), "File list must match the length of the primer list."
    # Defining the base directory
    file = file_list[0]  # Making the base directory that of the first file
    if platform == "linux" or platform == "linux2" or platform == "darwin":
        base_dir = "/".join(file.split("/")[:-1])
    elif platform == "win32":
        base_dir = "\\".join(file.split("\\")[:-1])
    else:
        base_dir = "/".join(file.split("/")[:-1])
    # Generating the new file for collating gene lists
    now = datetime.now()
    date = now.strftime("%y_%m_%d")
    file_name = f"collated_pooled_wish_list_{date}.txt"
    file_name = get_unique_filename(base_dir, file_name)
    file = f"{base_dir}/{file_name}"
    with open(file, 'w') as new_file:
        lib_num = 0
        for gene_file, primer_tuple in zip(file_list, primer_list):
            gene_list = read_gene_list(gene_file)
            left_primer, right_primer = primer_tuple
            for gene in gene_list:
                new_file.write(f"{gene}\t{left_primer}\t{right_primer}\t{lib_num}\n")
            lib_num += 1


def get_primer_list(primer_file: str,
                    num_sets: int) -> List[Tuple[str, str]]:
    """
    This function generates a list of primers

    Parameters
    ----------
    primer_file : str
        The primer file containing the primer lists.
    num_sets : int
        The number of primer pairs to pull.

    Returns
    -------
    primer_list: List[Tuple[str, str]]
        The list of primer tuples to use.
    """
    df = pd.read_csv(primer_file, sep="\t")
    # df = df.drop_duplicates(df)  # Making sure there are no duplicate primer pairs
    assert num_sets <= len(df), "Too many libraries are being collated for the number of primer sets available."
    rand_indices = np.random.choice(np.arange(0, len(df)), num_sets, replace=False)
    df = df.iloc[rand_indices]
    primer_list = list(zip(df['left'].values, df['right'].values))
    return primer_list


def main(raw_args=None):
    # Setting up CLI
    parser = ArgumentParser()
    parser.add_argument("--wishlist_files",
                        type=str,
                        nargs="+",
                        help="Path to txt file with gene wishlist.")
    parser.add_argument("--primer_file",
                        help="File that contains user-defined primers. OR None for randomized primers.",
                        default=None)
    # Parsing the Args
    args = parser.parse_args(raw_args)
    file_list = args.wishlist_files
    primer_file_ud = args.primer_file  # user-defined primer file

    # Getting primer file path
    file_path = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    file_path = os.path.join(file_path, '..', 'data')
    # Reading in list of default primers
    primer_file = os.path.join(file_path, 'pooled_primers.txt')
    # Defining the number of primer pairs to use
    num_pairs = len(file_list)
    # Checking if user-defined file with primers is given
    if primer_file_ud is None:
        primer_list = get_primer_list(primer_file, num_pairs)
    elif type(primer_file_ud) == str:
        primer_list = get_primer_list(primer_file_ud, num_pairs)
    else:
        raise Exception("The user provided file path for the primer file must have type str.")
    generate_pooled_list(file_list, primer_list)


if __name__ == "__main__":
    main()
