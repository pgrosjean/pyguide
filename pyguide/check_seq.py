from typing import List, Dict, Union
import numpy as np
import os
from pyguide.guide import open_txt_database, get_current_date, get_unique_filename
from argparse import ArgumentParser


def check_seq(file_list: List[str], organism: str, ai: str) -> List[str]:
    """
    This function takes the list of .seq files and checks them
    against all guides.

    Parameters
    ----------
    file_list: List[str]
        List of files containing sequencing reads.
    organism: str
        Whether to check against mouse or human guides.
    ai: str
        Whether CRISPRi or CRISPRa guides.

    Returns
    -------
    guide_id_list: np.ndarray
        Array of guide ids corresponding to the input sequences.
    """
    guide_map = get_map(organism=organism.lower(), ai_status=ai.lower())
    vec_read_seq = np.vectorize(read_seq_file)
    vec_get_protospacer = np.vectorize(find_protospacer)
    vec_get_guide_id = np.vectorize(guide_map.get)
    # vectorized operations
    seq_arr = vec_read_seq(np.array(file_list))
    protospacer_arr = vec_get_protospacer(seq_arr)
    guide_id_arr = vec_get_guide_id(protospacer_arr)
    return guide_id_arr


def read_seq_file(file: str) -> str:
    """
    This function reads in a .seq file.

    Parameters
    ----------
    file: str
        The seq file name.

    Returns
    -------
    seq: str
        The nucleotide sequence.
    """
    with open(file, 'r') as f:
        return "".join(str(f.read()).split("\n")[1:])


def get_map(organism: str, ai_status: str) -> Dict[str, str]:
    """
    This function generates a dictionary from

    Parameters
    ----------
    organism: str
        organism that guides were ordered for.
    ai_status: str
        Whether CRISPRa/i

    Returns
    -------
    guide_map: Dict[str, str]
        Map from protospacer to guide_id.
    """
    file_path = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    file_path = os.path.join(file_path, '..', 'data')
    if organism == "mouse":
        if ai_status == "i":
            file = os.path.join(file_path, "CRISPRi_mouse_v2.txt")
        elif ai_status == 'a':
            file = os.path.join(file_path, "CRISPRa_mouse_v2.txt")
        else:
            raise Exception("ai status must be either a or i")
    else:
        if ai_status == "i":
            file = os.path.join(file_path, "CRISPRi_v2.txt")
        elif ai_status == 'a':
            file = os.path.join(file_path, "CRISPRa_v2.txt")
        else:
            raise Exception("ai status must be either a or i")
    db_df = open_txt_database(file, header=False)  # Database DataFrame
    # Adding column names to database
    db_df = db_df.rename(columns={0: 'name', 1: 'gene', 2: 'score', 3: 'seq'})
    guide_map = dict(zip(db_df['seq'], db_df['name']))
    return guide_map


def find_protospacer(seq: str) -> Union[str, None]:
    """
    This function finds and extracts the protospacer from the seq file.

    Parameters
    ----------
    seq: str
        Sequence from seq file.

    Returns
    -------
    guide_id: Union[str, None]
        Guide id if found otherwise None.
    """
    # For pMK1334 and pLG15
    idx = seq.find("GTTTAAGAGCTAAGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTC")
    if idx == -1:
        return None
    else:
        protospacer_idx = idx - 20
        protospacer = seq[protospacer_idx:protospacer_idx + 20]
        return protospacer


def collate_files(directory: str, suffix: str) -> List[str]:
    """
    This function collates file with a user defined suffix.

    Parameters
    ----------
    directory: str
        Directory to walk recursively for files that end with suffix.
    suffix: str
        The suffix for filtering files.

    Returns
    -------
    file_list: List[str]
        List of files in directory that end with suffix.
    """
    file_list = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(suffix):
                file_list.append(os.path.join(root, file))
    return file_list


def main():
    parser = ArgumentParser()
    parser.add_argument("--file_dir",
                        type=str,
                        help="Path to directory with .seq files to check.")
    parser.add_argument("--ai",
                        type=str,
                        default="i",
                        help="i for CRISPRi or a for CRISPRa.")
    parser.add_argument("--organism",
                        type=str,
                        default="human",
                        help="Guides that target the mouse or human genome.")
    args = parser.parse_args()
    # Generating file list
    file_list = collate_files(args.file_dir, ".seq")
    assert len(file_list) != 0, "No .seq files found in directory."
    # Getting guide IDs for perfect matches
    guide_id_arr = check_seq(file_list, args.organism, args.ai)
    # Writing log file
    date = get_current_date()
    file_dir = args.file_dir
    log_file_name = f"log_check_seq_{date}.txt"
    log_file_name = os.path.join(file_dir, get_unique_filename(args.file_dir, log_file_name))
    with open(log_file_name, "w") as f:
        for file, guide_id in sorted(list(zip(file_list, list(guide_id_arr))), key=lambda x: x[1]):
            f.write(f"{file} \t {guide_id}\n")


