import numpy as np
import pandas as pd
from typing import List, Tuple
from sys import platform
from argparse import ArgumentParser
from pyguide.guide import get_unique_filename, read_gene_list

####################
# Defining Functions
####################


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
    # Generating the new file for collating gene lists
    file_name = "collated_pool_wish_list.txt"
    file_name = get_unique_filename(base_dir, file_name)
    file = f"{base_dir}/{file_name}"
    with open(file) as new_file:
        for gene_file, primer_tuple in zip(file_list, primer_list):
            gene_list = read_gene_list(gene_file)
            left_primer, right_primer = primer_tuple
            for gene in gene_list:
                new_file.write(f"{gene} \t {left_primer} \t {right_primer} \n")


def get_primer_list(primer_file: str,
                    num_sets: int,
                    seed: int = 42) -> List[Tuple[str, str]]:
    """
    This function generates a list of primers

    Parameters
    ----------
    primer_file : str
        The primer file containing the primer lists.
    num_sets : int
        The number of primer pairs to pull.
    seed : int
        The seed for randomly

    Returns
    -------
    primer_list: List[Tuple[str, str]]
        The list of primer tuples to use.
    """
    np.random.seed(seed=seed)
    df = pd.read_csv(primer_file, sep="\t")
    rand_indices = np.random.random_integers(0, len(df), num_sets)
    df = df.iloc[rand_indices]
    primer_list = list(zip(df['left'].values, df['right'].values))
    return primer_list


def main():
    # Setting up CLI
    parser = ArgumentParser()
    # TODO: CHECK TO MAKE SURE THIS WORKS FOR ADDING ALL FILES TOGETHER
    parser.add_argument("--wishlist_files",
                        type=list,
                        help="Path to txt file with gene wishlist.")
    parser.add_argument("--primer_list",
                        type=List[Tuple[str][str]],
                        help="List of user defined primers if wanted.",
                        default=None)
    parser.add_argument("--seed",
                        type=int,
                        help="Seed for random primer sequences if user_defined primer list not provided.",
                        default=None)
    # Parsing the Args
    args = parser.parse_args()
    file_list = args.wishlist_files
    primer_list = args.primer_list

    # Getting primer file path
    file_path = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    file_path = os.path.join(file_path, '..', 'data')
    primer_file = os.path.join(file_path, 'pooled_primers.txt')
    if primer_list is None:
        num_pairs = len(file_path)
        if args.seed is None:
            seed = np.random.randn()
        primer_list = get_primer_list(primer_file,
                                      num_pairs,
                                      seed)
    generate_pooled_list(file_list, primer_list)

