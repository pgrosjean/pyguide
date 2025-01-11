import pandas as pd
import numpy as np
import os
from typing import List, Dict, Tuple, Union
from argparse import ArgumentParser
from datetime import datetime
from sys import platform
import mygene


## TODO: Add a step for the plasmid backbone
## TODO: Check to make sure that we have the correct plasmid context for different backbones
## TODO: When outputting the check seq log file add a column with the plasmid number assigned to that correct sequence

####################
# Defining Functions
####################

####################
# I/O Functionality
####################

def read_wishlist(file: str) -> Tuple[List[str], List[str]]:
    """
    Reads a wishlist file containing either guide IDs or gene names.

    Parameters
    ----------
    file : str
        File containing the wishlist.

    Returns
    -------
    guide_ids : List[str]
        List of guide IDs from the wishlist.
    gene_names : List[str]
        List of gene names from the wishlist.
    """
    assert file.endswith(".txt"), "Wishlist file must be a .txt file."
    guide_ids, gene_names = [], []
    with open(file) as f:
        for line in f:
            item = line.strip()
            if "_" in item:  # Assume guide IDs contain underscores
                guide_ids.append(item)
            elif item:  # Non-empty lines without underscores are gene names
                gene_names.append(item)
    return guide_ids, gene_names


def read_gene_list(file: str) -> List[str]:
    """
    This function reads in gene names from a text file and returns a numpy array of gene names.

    Parameters
    ----------
        file : str
            The text file containing the string.

    Returns
    -------
        gene_list : List[str]
            The list of wishlist genes.
    """
    assert file.endswith(".txt"), "Gene files must be .txt file"
    gene_list = []
    with open(file) as f:
        lines = f.readlines()
        for line in lines:
            if line.strip() != '':
                gene_list.append(line.strip())
    return gene_list


def read_gene_list_pooled(file: str) -> Tuple[List[str], List[str], List[str], List[str]]:
    """
    This file reads in a pooled library gene wishlist file.

    Parameters
    ----------
    file : str
        The collated pooled library gene wishlist file.

    Returns
    -------
    genes : List[str]
    left_primers : List[str]
    right_primers : List[str]
    lib_num : List[str]
    """
    # reading in the gene list information as a DataFrame
    df = pd.read_csv(file, sep="\t", header=None)
    assert len(df.columns) == 4, f"Must run pyguide-collate on gene wishlists before pooled ordering. {df}"
    genes = list(df[0].values)
    left_primers = list(df[1].values)
    right_primers = list(df[2].values)
    lib_num = list(df[3].values)
    return genes, left_primers, right_primers, lib_num


####################
# Checking Database Functionality
####################

def open_txt_database(file: str, header: bool = True) -> pd.DataFrame:
    """
    This function reads in a txt database into a pandas dataframe.

    Parameters
    ----------
    file : str
        The database file.
    header : bool
        Whether the database has columns.

    Returns
    -------
    data_df : pd.DataFrame
        The dataframe of the database.
    """
    # handling header
    if header:
        head_row = 0
    else:
        head_row = None
    # reading database into a dataframe
    data_df = pd.read_csv(file, sep='\t', header=head_row)
    return data_df


def filter_database(df: pd.DataFrame,
                    guide_ids: List[str],
                    gene_names: List[str]) -> pd.DataFrame:
    """
    Filters the database by guide IDs or gene names.

    Parameters
    ----------
    df : pd.DataFrame
        The guide database dataframe.
    guide_ids : List[str]
        List of guide IDs to filter by.
    gene_names : List[str]
        List of gene names to filter by.

    Returns
    -------
    filtered_df : pd.DataFrame
        Filtered dataframe.
    """
    if len(guide_ids) > 0:
        filtered_df_guide = df[df['name'].isin(guide_ids)]
        assert len(filtered_df_guide) > 0, "No guides found in database."
    if len(gene_names) > 0:
        filtered_df_gene = df[df['gene'].isin(gene_names)]
        assert len(filtered_df_gene) > 0, "No genes found in database."
    if len(guide_ids) > 0 and len(gene_names) > 0:
        filtered_df = pd.concat((filtered_df_guide, filtered_df_gene), axis=0)
    elif len(guide_ids) > 0:
        filtered_df = filtered_df_guide
    elif len(gene_names) > 0:
        filtered_df = filtered_df_gene
    assert len(filtered_df) > 0, "No guides or genes found in database."
    return filtered_df


def get_guide_db(guide_ids: List[str],
                 gene_names: List[str],
                 ai_status: str,
                 organism: str) -> pd.DataFrame:
    """
    This function reads in the human guide database and filters for genes of interest from the user defined gene list.

    Parameters
    ----------
    guide_ids : List[str]
        The list of guide IDs.
    gene_list : List[str]
        The list of genes of interest.
    ai_status : str
        Whether to use CRISPRi/a.
    organism : str
        The name of the organism for generating guides.

    Returns
    -------
    db_df : pd.DataFrame
        The filtered dataframe of the database.
    """
    assert organism in ["mouse", "human"], "Organism must be either mouse or human."
    # Opening CRISPRi/a sgRNA database
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
            if len(gene_names) == 1 and gene_names[0] == 'negative_control':
                file = os.path.join(file_path, "CRISPRi_v2.txt")
            else:
                file = os.path.join(file_path, "CRISPRa_v2.txt")
        else:
            raise Exception("ai status must be either a or i")
    db_df = open_txt_database(file, header=False)  # Database DataFrame
    # Adding column names to database
    db_df = db_df.rename(columns={0: 'name', 1: 'gene', 2: 'score', 3: 'seq'})
    # Filtering full database for genes of interest
    db_df = filter_database(db_df, guide_ids, gene_names)
    return db_df


def get_empirical_db() -> pd.DataFrame:
    """
    This function reads in the empirical guide database.

    Returns
    -------
    empirical_df: pd.DataFrame
        Dataframe containing empirical phenotypes for guides in the databases.
    """
    file_path = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    file_path = os.path.join(file_path, '..', 'data')
    file = os.path.join(file_path, "sgRNAempiricalPhenotypes.csv")
    empirical_df = pd.read_csv(file)
    return empirical_df


def get_local_db(gene_list: List[str]) -> pd.DataFrame:
    """
    This function filters and returns the local Kampmann Lab Database for cloned sgRNAs.

    Parameters
    ----------
    gene_list : List[str]

    Returns
    -------
    kl_df : pd.DataFrame
        The local Kampmann Lab Database
    """
    file_path = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    file_path = os.path.join(file_path, '..', 'data')
    local_db_file = os.path.join(file_path, "human_sgrnas.txt")  # Local Kampmann Lab Database
    kl_df = open_txt_database(local_db_file, header=True)
    # Adding gene column to local Kampmann lab database
    kl_df["gene"] = kl_df["Short Name"].apply(lambda x: x.split("_")[0])
    # Filtering local database for genes of interest
    kl_df = filter_database(kl_df, [], gene_list)
    return kl_df


def get_cloned_status(local_df: pd.DataFrame,
                      guides_per_gene: int,
                      ai_status: str) -> Dict[str, int]:
    """
    This function determines whether guides have already been ordered for the genes in the user defined gene list and
    if more need to be ordered to reach the user defined guides_per_gene number.

    Parameters
    ----------
    local_df : pd.DataFrame
        Local Kampmann Lab database
    guides_per_gene : int
        Number of guides per gene to order
    ai_status : str
        Whether CRISPRi or CRISPRa

    Returns
    -------
    cloned: Dict[str, int]
        Dictionary containing the number of guides that have already been cloned be cloned for each gene that has
        already had a guide cloned.
    """
    cloned = {}
    if local_df.shape[0] == 0:
        return cloned
    else:
        # Filtering to only apply to either CRISPRa or CRISPRi
        local_df['ai_status'] = local_df['Short Name'].apply(lambda x: x.split('_')[-1][-2])
        local_df = local_df.loc[local_df['ai_status'] == ai_status]
        for gene_name, gene_df in local_df.groupby('gene'):
            if gene_df.shape[0] >= guides_per_gene:
                cloned[gene_name] = guides_per_gene
            else:
                cloned[gene_name] = gene_df.shape[0]
    return cloned


def collate_guides(guides_per_gene: int,
                   hg_db: pd.DataFrame) -> pd.DataFrame:
    """
    This function collates guides from the human sgRNA Database. The hg_db is either for CRISPRi or CRISPRa.

    Parameters
    ----------
    guides_per_gene : int
        Number of guides per gene.
    hg_db : pd.DataFrame
        Human sgRNA database data frame.

    Returns
    -------
    collated_df : pd.DataFrame
        Data frame with collated guides.
    """
    # Taking the top guides_per_gene guides for ordering
    collated_df = hg_db.sort_values('score', ascending=False).groupby('gene').head(guides_per_gene)
    collated_df = collated_df.sort_values('gene')
    collated_df = collated_df.reset_index(drop=True)
    return collated_df


def groupby_and_rank(df, col_1, col_2, col_3, head):
    """
    Group a dataframe by col_3, then rank the resulting groups by col_1 (if
    possible) or col_2 (if not). The head parameter determines how many rows
    from each group are returned.

    Parameters
    ----------
    df: pd.DataFrame
        Input dataframe.
    col_1: str
        Empirical guide score column name.
    col_2: str
        Computational guide score column name.
    col_3: str
        Group by column name.
    head: int
        number of rows to take from each group.

    Returns
    -------
    pd.DataFrame
        Ranked and filtered dataframe.
    """

    def rank_col_1(group):
        """
        Helper function for ranking based on col_1. Assumes that the group
        contains both numerical values and NaNs in col_1, and numerical
        values in col_2.
        """
        col_1_mask = group[col_1].notna()
        num_col_1 = col_1_mask.sum()

        col_2_ranks = group.loc[~col_1_mask, col_2].apply(lambda x: -x).argsort().argsort()
        if num_col_1 == 0:
            return col_2_ranks
        col_1_ranks = np.empty(len(group))
        col_1_ranks.fill(np.nan)
        col_1_ranks[col_1_mask] = (-group.loc[col_1_mask, col_1]).argsort().argsort()
        col_1_ranks[np.isnan(col_1_ranks)] = col_2_ranks + num_col_1
        assert np.sum(np.isnan(col_1_ranks)) == 0
        return col_1_ranks

    groups = df.groupby(col_3)

    ranked = []
    for _, group in groups:
        ranks = rank_col_1(group)
        ranked_group = group.assign(rank=ranks)
        ranked_group = ranked_group.sort_values(by="rank", na_position="last")
        ranked_group = ranked_group.head(head)
        ranked.append(ranked_group)

    result = pd.concat(ranked, axis=0)
    result = result.drop(columns=["rank"])
    result = result.reset_index(drop=True)
    assert result["score"].isna().sum() == 0
    return result


def collate_guides_empirical(guides_per_gene: int,
                             db: pd.DataFrame,
                             screen_db: pd.DataFrame) -> pd.DataFrame:
    """
    This function collates guides from the provided sgRNA Database using the empirical
    scores first and then using the computationally predicted scores second.

    Parameters
    ----------
    guides_per_gene: int
        The number of guides to order per gene.
    db: pd.DataFrame
        The sgRNA database.
    screen_db: pd.DataFrame
        The database of empirical guide RNA phenotype scores.

    Returns
    -------
    pd.DataFrame
        The filtered dataframe with the specific guides for ordering.
    """
    assert db["score"].isna().sum() == 0
    # Merge the two dataframes
    df = db.merge(screen_db, left_on='name', right_on='sgRNA', how='left')
    # Call the groupby_and_rank function
    filtered_df = groupby_and_rank(df, 'meanEmpiricalPhenotype', 'score', 'gene', guides_per_gene)
    collated_df = filtered_df.loc[:, db.columns]
    collated_df = collated_df.sort_values('gene')
    collated_df = collated_df.reset_index(drop=True)
    # return the filtered dataframe
    return collated_df


def collate_ntcs(num_ntc: int, db: pd.DataFrame):
    """
    This function collates a DataFrame with a set of num_ntc number of randomly selected NTCs.

    Parameters
    ----------
    num_ntc: int
        The number of NTC guides to order
    db: pd.DataFrame
        The sgRNA data frame.

    Returns
    -------
    collated_df: pd.DataFrame
        The collated dataframe with the non-targeting control guides.
    """
    if num_ntc <= 2:
        num_ntc = 2
    assert num_ntc >= 0
    db_df = db.loc[db['gene'] == 'negative_control']
    collated_df = db_df.sample(n=num_ntc)
    return collated_df


def filter_cloned_guides(collated_df: pd.DataFrame,
                         cloned_dict: Dict[str, int]) -> pd.DataFrame:
    """
    This function filters out any of the previously cloned guides so that only non-cloned guides are ordered.

    Parameters
    ----------
    collated_df : pd.DataFrame
        The collated df from the sgRNA Data base
    cloned_dict : Dict[str, int]

    Returns
    -------
    filtered_df : pd.DataFrame
        The filtered data frame.
    """
    filtered_df = collated_df.copy()
    filtered_df = filtered_df.sort_values('score', ascending=False).sort_values('gene')
    # Setting up boolean indexing vector
    bool_ind = np.ones(filtered_df.shape[0])
    # Changing boolean indexing vector to account for already cloned guides
    if len(filtered_df) == 0:
        filtered_df = filtered_df.reset_index(drop=True)
        return filtered_df
    else:
        for gene_name, num in cloned_dict.items():
            gene_filter = filtered_df[filtered_df['gene'] == gene_name].sort_values('score', ascending=False)
            gene_cond = gene_filter.head(num).index.values
            bool_ind[gene_cond] = 0
        bool_ind = bool_ind.astype(bool)
        filtered_df = filtered_df.iloc[bool_ind]
        filtered_df = filtered_df.reset_index(drop=True)
        return filtered_df


####################
# Helper Functions
####################

def get_current_date() -> str:
    """
    This function returns the current date.

    Returns
    -------
    date : str
        The current date.
    """
    now = datetime.now()
    date = now.strftime("%y_%m_%d")
    return date


def reverse_compliment(seq: str) -> str:
    """
    This function returns the reverse compliment of a DNA sequence.

    Parameters
    ----------
    seq : str
        The DNA sequence for reverse compliment.

    Returns
    -------
    rev_comp : str
        reverse compliment of input DNA sequence.
    """
    rev_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    rev_seq = seq[::-1]
    rev_comp = ""
    for nuc in rev_seq:
        rev_comp += rev_dict[nuc]
    return rev_comp


def get_short_names(df: pd.DataFrame,
                    ai_status: str) -> List[str]:
    """
    This function generates the short names for guides. e.g. STAT3_i1 for the top ranked CRISPRi guide.

    Parameters
    ----------
    df : pd.DataFrame
        Data frame for finding short names.
    ai_status : str
        Whether CRISPRi or CRISPRa.

    Returns
    -------
    short_names : List[str]
    """
    assert df['score'].isna().sum() == 0
    guide_ranks = df.groupby('gene')['score'].rank(ascending=False).astype('int').astype(str)
    short_names = list(df['gene'] + '_' + ai_status + guide_ranks)
    return short_names


def split_dataframe(df: pd.DataFrame,
                    chunk_size: int = 96) -> List[pd.DataFrame]:
    """
    This function splits a dataframe into chunks of user defined size.

    Parameters
    ----------
    df : pd.DataFrame
        Data frame to chunk.
    chunk_size : int
        size of chunks to make.

    Returns
    -------
    chunks : List[pd.DataFrame]
        Chunked data frame
    """
    chunks = list()
    num_chunks = len(df) // chunk_size + 1
    for i in range(num_chunks):
        chunks.append(df[i * chunk_size:(i + 1) * chunk_size])
    return chunks


def get_unique_filename(base_dir: str,
                        filename: str) -> str:
    """
    This function checks if the filename already exists at a directory and if it does it returns and unique updated
    file name.

    Parameters
    ----------
    base_dir : str
        Directory to check for file uniqueness.
    filename : str
        The file name to check.

    Returns
    -------
    unique_filename : str
        The unique filename for the directory.
    """
    file_list = os.listdir(base_dir)
    c = 1
    fname_check = filename
    while fname_check in file_list:
        split_file = filename.split('.')
        assert len(split_file) == 2, "Do not use a . in your filenames."
        fname = split_file[0]
        suffix = split_file[-1]
        fname_check = fname + f"_{c}." + suffix
        c += 1
    unique_filename = fname_check
    return unique_filename


####################
# Arrayed Ordering Functions
####################

def write_arrayed_csv(collated_df: pd.DataFrame,
                      ai_status: str,
                      name: str,
                      base_dir: str):
    """
    This function writes the ordering csv file used to order guides for arrayed guide cloning.

    Parameters
    ----------
    collated_df : pd.DataFrame
        The dataframe with the collated guides that need to be ordered.
    ai_status : str
        Whether CRISPRi or CRISPRa.
    name : str
        Name of user.
    base_dir : str
        Directory to save csv file to.

    Returns
    -------
    None
    """
    # Resetting the index of the collated data frame
    collated_df = collated_df.reset_index(drop=True)

    # Chunking collated data frame for writing to DNA plate csv files
    chunked_df = split_dataframe(collated_df, chunk_size=96)

    # Defining well ids
    well_ids = []
    for letter in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:
        for num in np.arange(1, 13):
            num = str(num)
            well_ids.append(letter + num)

    # Iterating through each plate and writing csv file
    for plate_num, plate_df in enumerate(chunked_df):
        plate_num = plate_num + 1
        bp_name = []  # bottom plate names
        bp_oligos = []  # bottom plate oligos
        tp_name = []  # top plate names
        tp_oligos = []  # top plate oligos
        pw_ids = []  # this plates well ids
        short_names = get_short_names(plate_df, ai_status)
        for guide_name, proto_spacer, well_id in zip(short_names, plate_df['seq'], well_ids):
            # defining protospacer
            proto_spacer = proto_spacer.upper()
            # generating oligonucleotide sequences
            top_oligo = "TTG" + proto_spacer + "GTTTAAGAGC"
            bottom_oligo = "TTAGCTCTTAAAC" + reverse_compliment(proto_spacer) + "CAACAAG"
            # appending to corresponding plates
            tp_name.append(guide_name + "_top")
            bp_name.append(guide_name + "_bottom")
            bp_oligos.append(bottom_oligo)
            tp_oligos.append(top_oligo)
            pw_ids.append(well_id)
        # Generating DataFrames to write to csv files
        tp_df = pd.DataFrame({'Well Position': pw_ids,
                              'Name': tp_name,
                              'Sequence': tp_oligos})
        bp_df = pd.DataFrame({'Well Position': pw_ids,
                              'Name': bp_name,
                              'Sequence': bp_oligos})
        # Writing dataframes to csv files
        date = get_current_date()
        tp_filename = get_unique_filename(base_dir, f"order_arrayed_{name}_{date}_top_plate_{plate_num}.csv")
        bp_filename = get_unique_filename(base_dir, f"order_arrayed_{name}_{date}_bottom_plate_{plate_num}.csv")
        if platform == "linux" or platform == "linux2" or platform == "darwin":
            tp_df.to_csv(f"{base_dir}/{tp_filename}",
                         index=False)
            bp_df.to_csv(f"{base_dir}/{bp_filename}",
                         index=False)
        elif platform == "win32":
            tp_df.to_csv(f"{base_dir}\\{tp_filename}",
                         index=False)
            bp_df.to_csv(f"{base_dir}\\{bp_filename}",
                         index=False)


def write_basic_log_file(missing_genes: List[str],
                         changed_genes: Dict[str, str],
                         name: str,
                         base_dir: str,
                         order_format: str):
    """
    This file writes a basic log file.

    Parameters
    ----------
    missing_genes : List[str]
        Any genes in wish list that were not found in the overall database.
    changed_genes : None or Dict[str, str]
        A dictionary containing changes from queried genes to actual gene name in database or None.
    name : str
        Name of user.
    base_dir : str
        The directory to save the log file to.
    order_format : str
        The order format of either arrayed or single.

    Returns
    -------
    None
    """
    # Setting up Log File
    date = get_current_date()
    filename = get_unique_filename(base_dir, f"log_file_{order_format}_{name}_{date}.txt")
    if platform == "linux" or platform == "linux2" or platform == "darwin":
        file = f"{base_dir}/{filename}"
    elif platform == "win32":
        file = f"{base_dir}\\{filename}"
    else:
        file = f"{base_dir}/{filename}"

    # Writing to Log file
    with open(file, 'w') as f:
        f.write(f"This query was performed by {name} on {date} \n \n")
        if len(missing_genes) != 0:
            f.write("\n")
            for missing_gene in missing_genes:
                f.write(f"Gene {missing_gene} was not found in the DataBase." + "\n")
            f.write("\n")
        if changed_genes is not None:
            for changed_gene in changed_genes.keys():
                f.write(f"Query for {changed_gene} was mapped to {changed_genes[changed_gene]} in database." + "\n")
        else:
            f.write("All queried guides have not yet been cloned and no errors were encountered." + "\n")


####################
# Single (one-at-a-time) Ordering Functions
####################

def write_single_csv(filtered_df: pd.DataFrame,
                     ai_status: str,
                     name: str,
                     base_dir: str):
    """
    This function writes the ordering csv file used to order guides for single guide cloning.

    Parameters
    ----------
    filtered_df : pd.DataFrame
        The dataframe with the collated guides that need
        to be ordered.
    ai_status : str
        Whether CRISPRi or CRISPRa.
    name : str
        Name of user.
    base_dir : str
        Directory to save csv file to.

    Returns
    -------
    None
    """
    short_names = get_short_names(filtered_df, ai_status)
    name_list = []
    oligo_list = []
    # Writing top and bottom oligonucleotides
    for guide_name, proto_spacer in zip(short_names, filtered_df['seq']):
        name_list.append(guide_name + '_top')
        name_list.append(guide_name + '_bottom')
        proto_spacer = proto_spacer.upper()
        top_oligo = "TTG" + proto_spacer + "GTTTAAGAGC"
        bottom_oligo = "TTAGCTCTTAAAC" + reverse_compliment(proto_spacer) + "CAACAAG"
        oligo_list.append(top_oligo)
        oligo_list.append(bottom_oligo)
    scale = ['25nm'] * len(oligo_list)
    purification = ['STD'] * len(oligo_list)
    order_df = pd.DataFrame({'Name': name_list,
                             'Sequence': oligo_list,
                             'Scale': scale,
                             'Purification': purification})
    date = get_current_date()
    filename = get_unique_filename(base_dir, f"order_single_{name}_{date}.csv")
    if platform == "linux" or platform == "linux2" or platform == "darwin":
        order_df.to_csv(f"{base_dir}/{filename}", index=False)
    elif platform == "win32":
        order_df.to_csv(f"{base_dir}\\{filename}", index=False)


def write_single_log_file(local_df: pd.DataFrame,
                          missing_genes: List[str],
                          changed_genes: Dict[str, str],
                          name: str,
                          ai_status: str,
                          base_dir: str):
    """
    Writes log file of where to find previously cloned guides and any error logging.

    Parameters
    ----------
    local_df : pd.DataFrame
        The dataframe of local library database.
    changed_genes : None or Dict[str, str]
        A dictionary containing changes from queried genes to actual gene name in database or None.
    missing_genes : List[str]
        Any genes in wish list that were not found in the overall database.
    name : str
        Name of user.
    ai_status : str
        CRISPRa or CRISPRi.
    base_dir : str
        The directory to save the log file to.

    Returns
    -------
    None
    """
    # Filtering for CRIPSRi/a status first
    local_df['ai_status'] = local_df['Short Name'].apply(lambda x: x.split('_')[-1][0])
    local_df = local_df.loc[local_df['ai_status'] == ai_status]

    # Setting up Log File
    date = get_current_date()
    filename = get_unique_filename(base_dir, f"log_file_{name}_{date}.txt")
    if platform == "linux" or platform == "linux2" or platform == "darwin":
        file = f"{base_dir}/{filename}"
    elif platform == "win32":
        file = f"{base_dir}\\{filename}"
    else:
        file = f"{base_dir}/{filename}"  # always default to linux and macos

    # Writing to Log file
    with open(file, 'w') as f:
        f.write(f"This query was performed by {name} on {date} \n \n")
        if local_df.shape[0] != 0:
            local_df = local_df.sort_values('Short Name')
            for short_name, number in zip(local_df['Short Name'].values, local_df['#'].values):
                f.write(f"{short_name} has been cloned with the ID: {number}" + "\n")
        if len(missing_genes) != 0:
            f.write("\n")
            for missing_gene in missing_genes:
                f.write(f"Gene {missing_gene} was not found in the DataBase." + "\n")
            f.write("\n")
        if changed_genes is not None:
            for changed_gene in changed_genes.keys():
                f.write(f"Query for {changed_gene} was mapped to {changed_genes[changed_gene]} in database." + "\n")
        if local_df.shape[0] == 0 and len(missing_genes) == 0:
            f.write("All queried guides have not yet been cloned and no errors were encountered." + "\n")


####################
# Pooled Ordering Functions
####################

def write_pooled_txt(collated_df: pd.DataFrame,
                     name: str,
                     base_dir: str,
                     batch_retest: bool = False):
    """
    This function writes a text file containing the information for ordering pooled guide libraries from Agilent.

    Parameters
    ----------
    collated_df : pd.DataFrame
        Data frame with collated guides.
    name : str
        Name of user.
    base_dir : str
        The base directory to save text file to.

    Returns
    -------
    None
    """
    if batch_retest:
        date = get_current_date()
        filename = get_unique_filename(base_dir, f"order_batch_retest_{name}_{date}.txt")
        filename_2 = get_unique_filename(base_dir, f"order_batch_retest_{name}_{date}_info.csv")
    else:
        date = get_current_date()
        filename = get_unique_filename(base_dir, f"order_pooled_{name}_{date}.txt")
        filename_2 = get_unique_filename(base_dir, f"order_pooled_{name}_{date}_info.csv")
    if platform == "linux" or platform == "linux2" or platform == "darwin":
        file = f"{base_dir}/{filename}"
        file_2 = f"{base_dir}/{filename_2}"
    elif platform == "win32":
        file = f"{base_dir}\\{filename}"
        file_2 = f"{base_dir}\\{filename_2}"
    else:
        file = f"{base_dir}/{filename}"  # Always baseline assumes Linux OS
        file_2 = f"{base_dir}/{filename_2}"
    res_left = "CCACCTTGTTG"  # BstXI restriction enzyme
    res_right = "GTTTAAGAGCTAAGCTGG"  # Bpi1102I restriction enzyme
    with open(file, 'w') as f:
        seqs = collated_df['seq'].values
        left_primers = collated_df['left_primers'].values
        right_primers = collated_df['right_primers'].values
        for seq, left_primer, right_primer in zip(seqs, left_primers, right_primers):
            left_adapter = left_primer
            right_adapter = reverse_compliment(right_primer)
            # 5' --> 3'
            full_seq = left_adapter + res_left + seq + res_right + right_adapter
            full_seq_rev_comp = reverse_compliment(full_seq)
            f.write(full_seq + "\n")
            f.write(full_seq_rev_comp + "\n")
    collated_df.to_csv(file_2)


def write_pooled_log_file(missing_genes: List[str],
                          changed_genes: Dict[str, str],
                          name: str,
                          base_dir: str,
                          primer_df: pd.DataFrame):
    """
    This function writes a log file for pooled ordering that contains information about if a requested gene could not
    be found. If the requested gene was an alias or had a different name than the one that was originally queried this
    file contains that information. And additionally this file gives information about the primer sequences used for
    the different requested libraries.

    Parameters
    ----------
    missing_genes : List[str]
        Any genes in wish list that were not found in the overall database.
    changed_genes : None or Dict[str, str]
        A dictionary containing changes from queried genes to actual gene name in database or None.
    name : str
        Name of user.
    base_dir : str
        The directory to save the log file to.
    primer_df : pd.DataFrame
        The data frame containing the information about primers and library numbers.

    Returns
    -------
    None

    """
    # Setting up Log File
    date = get_current_date()
    filename = get_unique_filename(base_dir, f"log_file_pooled_{name}_{date}.txt")
    if platform == "linux" or platform == "linux2" or platform == "darwin":
        file = f"{base_dir}/{filename}"
    elif platform == "win32":
        file = f"{base_dir}\\{filename}"
    else:
        file = f"{base_dir}/{filename}"

    # Filtering primer_df for unique library and primer information
    primer_df = primer_df.drop(columns='gene_symbol')
    primer_df = primer_df.drop_duplicates()
    left_primers = list(primer_df['left_primers'].values)
    right_primers = list(primer_df['right_primers'].values)
    library_numbers = list(primer_df['lib_num'].values)

    # Writing to Log file
    with open(file, 'w') as f:
        f.write(f"This query was performed by {name} on {date} \n \n")
        if len(missing_genes) != 0:
            f.write("\n")
            for missing_gene in missing_genes:
                f.write(f"Gene {missing_gene} was not found in the DataBase." + "\n")
            f.write("\n")
        if changed_genes is not None:
            for changed_gene in changed_genes.keys():
                f.write(f"Query for {changed_gene} was mapped to {changed_genes[changed_gene]} in database." + "\n")
        else:
            f.write("All queried guides have not yet been cloned and no errors were encountered." + "\n")
        f.write(f"\nPrimer Information\n------------------\nLibrary Number\t Primer 1\tPrimer 2\n")
        for left_primer, right_primer, library_number in zip(left_primers, right_primers, library_numbers):
            f.write(f"{library_number}\t{left_primer}\t{right_primer}\n")


def write_batch_retest_log_file(name, base_dir, primer_df):
    """
    This function is used to write a log file for batch retesting of pooled libraries.

    Parameters
    ----------
    name : str
        Name of user.
    base_dir : str
        The directory to save the log file to.
    primer_df : pd.DataFrame
        The data frame containing the information about primers and library numbers.
    """
    # Setting up Log File
    date = get_current_date()
    filename = get_unique_filename(base_dir, f"log_file_batch_retest_{name}_{date}.txt")
    if platform == "linux" or platform == "linux2" or platform == "darwin":
        file = f"{base_dir}/{filename}"
    elif platform == "win32":
        file = f"{base_dir}\\{filename}"
    else:
        file = f"{base_dir}/{filename}"
    primer_df = primer_df.drop(columns='guide_id')
    primer_df = primer_df.drop_duplicates()
    left_primers = list(primer_df['left_primers'].values)
    right_primers = list(primer_df['right_primers'].values)
    library_numbers = list(primer_df['lib_num'].values)

    with open(file, 'w') as f:
        f.write(f"This batch retest library was generated by {name} on {date} \n \n")
        f.write(f"\nPrimer Information\n------------------\nLibrary Number\t Primer 1\tPrimer 2\n")
        for left_primer, right_primer, library_number in zip(left_primers, right_primers, library_numbers):
            f.write(f"{library_number}\t{left_primer}\t{right_primer}\n")


####################
# MyGene Querying and Symbol Mapping Functions
####################

def make_query_map(query_df: pd.DataFrame) -> Dict[str, str]:
    """
    This function generates a map between the queried name and the database name.

    Parameters
    ----------
    query_df: pd.DataFrame
        DataFrame containing information about the queried genes vs those that are in the database.

    Returns
    -------
    query_map: Dict[str, str]
    """
    query_map = {}
    for query, symbol, alias in zip(query_df.index, query_df['symbol'].values, query_df['alias'].values):
        if type(alias) != float:
            if type(alias) != list:
                alias = list(alias)
            query_map[symbol] = query
            for al in alias:
                query_map[al] = query
    return query_map


def query_genes(symbol_list: List[str], species: str) -> Union[pd.DataFrame, None]:
    """
    This function queries human genes using the mygenes API.

    Parameters
    ----------
    species : str
        The name of the species that you are querying (must be mouse or human)
    symbol_list : List[str]
        The list of gene symbols to match to other aliases and ensemble gene IDs and entrez IDs.

    Returns
    -------
    gene_df : Union[pd.DataFrame, None]
        Dataframe containing all the gene information, unless the gene name queried is not found where None is returned.
    """
    species = species.lower()
    if species == "human":
        species_num = 9606
    elif species == "mouse":
        species_num = 10090
    else:
        species_num = None
    assert species_num is not None, "Need to pass either human or mouse for species!"
    mg = mygene.MyGeneInfo()
    gene_map = mg.querymany(symbol_list,
                            scopes=['symbol', 'alias'],
                            species=species_num,
                            fields='ensembl.gene, entrezgene, symbol, alias',
                            as_dataframe=True,
                            returnall=True,
                            verbose=False)
    if len([x for x in gene_map['missing']]) == len(symbol_list):
        return None
    elif 'alias' not in gene_map['out'].columns:
        return None
    else:
        gene_df = gene_map['out']
        if 'notfound' in gene_df.columns:
            gene_df = gene_df.loc[gene_df['notfound'] != True, :]
        return gene_df


def get_guides_direct(guide_ids: List[str],
                      ai_status: str,
                      organism: str) -> pd.DataFrame:
    """
    This function gets guides directly from the database.

    Parameters
    ----------
    guide_ids : List[str]
        List of guide IDs to get from the database.
    ai_status : str
        Whether CRISPRi or CRISPRa.
    organism : str
        Whether to get guides that target the mouse or human.

    Returns
    -------
    guide_df : pd.DataFrame
        Data frame with the guides that were directly queried.
    """
    guide_df = get_guide_db(guide_ids, gene_names=[], ai_status=ai_status, organism=organism)
    return guide_df

####################
# Main entry point for all guide ordering
####################

def order_guides(guide_ids: List[str],
                 gene_names: List[str],
                 name: str,
                 ai_status: str,
                 guides_per_gene: int,
                 order_format: str,
                 base_dir: str,
                 check_db: bool,
                 organism: str,
                 ntc_frac: float = 0.2,
                 primer_df: Union[pd.DataFrame, None] = None):
    """
    Orders guides for cloning, supporting guide IDs and gene names.

    Parameters
    ----------
    guide_ids: List[str]
        List of guide IDs to order.
    gene_names: List[str]
        Names of genes for ordering.
    name: str
        Name of individual running script.
    ai_status: str
        Whether experiment is CRISPRa or CRISPRi.
    guides_per_gene: int
        Number of guides per gene.
    order_format: str
        Whether ordering arrayed, pooled, or single guides.
    base_dir: str
        The location to save the output files to.
    check_db: bool
        Whether the user is in the Kampmann Lab. Controls access to DataBase.
    organism: str
        Whether to get guides that target the mouse or human.
    primer_df: pd.DataFrame
        Data frame containing information about primers for pooled ordering.
    ntc_frac: float
        The number of ntc guides to order if pooled ordering.
    """
    # Making all str inputs for if else statements lower case
    ai_status = ai_status.lower()
    order_format = order_format.lower()
    possible_order_formats = ['arrayed', 'single', 'pooled', 'batch-retest']
    assert order_format in possible_order_formats, f"{order_format} not in {possible_order_formats}"
    assert guides_per_gene <= 10, "No more than 10 guides per gene allowed."

    # Pulling in sgRNA filtered database
    missing_genes = []
    mod_query_map = None
    if len(gene_names) > 0:
        db_df = get_guide_db(guide_ids=[], gene_names=gene_names, ai_status=ai_status, organism=organism)
        # Logging any genes in wishlist not found in database
        missing_genes = [gene for gene in gene_names if gene not in set(db_df['gene'].unique())]
        if len(missing_genes) > 0:
            # Checking missing genes against standardized web API for known genes
            missing_df = query_genes(missing_genes, species=organism)
            if missing_df is not None:
                query_map = make_query_map(missing_df)
                alias_names = np.hstack(list(query_map.keys()))
                symbols = np.hstack(missing_df['symbol'].values)
                missing_filter = list(np.hstack([alias_names, symbols]))
                missing_db_df = get_guide_db(guide_ids=[], gene_names=missing_filter, ai_status=ai_status, organism=organism)
                # ^ Returns any of aliases or symbols. So then I need to map these outputs to the original gene name
                # that was requested by the user
                db_df = pd.concat([missing_db_df, db_df])
                mapped_missed_genes = [query_map.get(gene) for gene in missing_db_df['gene'].values]
                mapped_missed_genes = [x for x in mapped_missed_genes if x is not None]
                all_mapped_genes = np.hstack([mapped_missed_genes, db_df['gene'].values])
                missing_genes = [gene for gene in gene_list if gene not in all_mapped_genes]
                # Capturing modified queries for returning to user
                mod_query_map = {}  # maps original query to updated name
                for gene_name in missing_db_df['gene']:
                    if gene_name in set(list(query_map.keys())):
                        mod_query_map[query_map[gene_name]] = gene_name
            else:
                mod_query_map = None
        else:
            mod_query_map = None
        # Updating gene_list
        gene_list = list(db_df['gene'].values)
        assert db_df['score'].isna().sum() == 0

    assert len(gene_names) > 0 or len(guide_ids) > 0, "Must provide either gene names or guide IDs for ordering."

    # Adding in functionality for ordering guides by names
    if len(gene_names) == 0:
        collated_df = get_guides_direct(guide_ids, ai_status=ai_status, organism=organism)
        collated_df.reset_index(drop=True)
    elif len(guide_ids) == 0:
        # Reading in empirical guide data
        empirical_df = get_empirical_db()
        # Collate guides
        collated_df = collate_guides_empirical(guides_per_gene, db_df, empirical_df)
        collated_df.reset_index(drop=True)
    else:
        guide_only_df = get_guides_direct(guide_ids, ai_status, organism)
        empirical_df = get_empirical_db()
        collated_df = collate_guides_empirical(guides_per_gene, db_df, empirical_df)
        collated_df = pd.concat([collated_df, guide_only_df], axis=0)
        collated_df.reset_index(drop=True)

    # Single Guide Ordering
    if order_format == "single":
        if check_db is True and organism == "human":
            # Pulling in local Lab Database
            local_df = get_local_db(gene_list)

            # Check if guides have been cloned before
            cloned_dict = get_cloned_status(local_df,
                                            guides_per_gene=guides_per_gene,
                                            ai_status=ai_status)

            # Filtering for previously cloned guides
            filtered_df = filter_cloned_guides(collated_df, cloned_dict)

            # Writing csv files for ordering guides
            write_single_csv(filtered_df,
                             ai_status=ai_status,
                             name=name,
                             base_dir=base_dir)

            # Writing log file for finding already cloned guides and logging errors
            write_single_log_file(local_df,
                                  missing_genes=missing_genes,
                                  changed_genes=mod_query_map,
                                  name=name,
                                  ai_status=ai_status,
                                  base_dir=base_dir)
        else:
            # Writing the csv file for ordering single guides with no filtering
            write_single_csv(collated_df,
                             ai_status=ai_status,
                             name=name,
                             base_dir=base_dir)

            # Writing the log file for the single guide ordering
            write_basic_log_file(missing_genes, mod_query_map, name, base_dir, "single")

    # Arrayed Guide Ordering
    elif order_format == "arrayed":
        # Writing the arrayed ordering csv file
        write_arrayed_csv(collated_df,
                          ai_status=ai_status,
                          name=name,
                          base_dir=base_dir)

        # Writing the arrayed log file
        write_basic_log_file(missing_genes=missing_genes,
                             changed_genes=mod_query_map,
                             name=name,
                             base_dir=base_dir,
                             order_format="arrayed")

    # Pooled guide ordering
    elif order_format == "pooled":
        assert primer_df is not None, "must provide primer information when pooled ordering."
        # changing gene names in the primer df if necessary
        if mod_query_map is not None:
            primer_df['gene_symbol'] = primer_df['gene_symbol'].apply(lambda x: mod_query_map.get(x, x))
        collated_df = collated_df.merge(primer_df, left_on='gene', right_on='gene_symbol')
        # Adding in negative controls to pooled libraries
        unique_library_info = primer_df.loc[:, ['left_primers', 'right_primers', 'lib_num']].drop_duplicates()
        unique_library_info = unique_library_info.reset_index(drop=True)
        ntc_db = get_guide_db([], ['negative_control'], ai_status=ai_status, organism=organism)
        ntc_sets = {}

        # Iterating through the libraries for ordering and generating NTCs
        for library_number in primer_df['lib_num'].unique():
            assert primer_df.loc[primer_df['lib_num'] == library_number].shape[0] >= 0
            num_ntc = int(ntc_frac * primer_df.loc[primer_df['lib_num'] == library_number].shape[0])
            if num_ntc <= 2:
                num_ntc = 2
            assert num_ntc >= 2, f"{num_ntc} Must have at least one non-targeting control when ordering pooled libraries."
            ntc_df = collate_ntcs(num_ntc=num_ntc, db=ntc_db)
            ntc_df['gene_symbol'] = ntc_df['gene']
            ntc_df = ntc_df.reset_index(drop=True)
            for_concat = np.repeat(unique_library_info.loc[unique_library_info['lib_num'] == library_number].values, len(ntc_df), axis=0)
            for_concat = pd.DataFrame({'left_primers': for_concat[:, 0],
                                       'right_primers': for_concat[:, 1],
                                       'lib_num': for_concat[:, 2]})
            ntc_df = pd.concat([ntc_df, for_concat], axis=1)
            ntc_df['lib_num'] = [library_number]*len(ntc_df)
            ntc_sets[library_number] = ntc_df

        for ntc_set in ntc_sets.values():
            collated_df = pd.concat([collated_df, ntc_set], axis=0)
        # Writing the pooled ordering txt file
        write_pooled_txt(collated_df, name, base_dir)
        write_pooled_log_file(missing_genes, mod_query_map, name, base_dir, primer_df)
    
    # Batch retest
    elif order_format == "batch-retest":
        collated_df = collated_df.merge(primer_df, left_on='name', right_on='guide_id')
        ntc_db = get_guide_db([], ['negative_control'], ai_status=ai_status, organism=organism)
        num_ntc = int(ntc_frac * len(primer_df))
        if num_ntc < 1:
            print(f"Not including any NTCs based on ntc_frac provided...")
        ntc_df = collate_ntcs(num_ntc=num_ntc, db=ntc_db)
        ntc_df = ntc_df.reset_index(drop=True)
        ntc_df['guide_id'] = ntc_df['gene']
        for_concat = pd.DataFrame({'left_primers': [primer_df['left_primers'].values[0]]*len(ntc_df),
                                      'right_primers': [primer_df['right_primers'].values[0]]*len(ntc_df),
                                      'lib_num': [primer_df['lib_num'].values[0]]*len(ntc_df)})
        ntc_df = pd.concat([ntc_df, for_concat], axis=1)
        collated_df = pd.concat([collated_df, ntc_df], axis=0)
        write_pooled_txt(collated_df, name, base_dir, batch_retest=True)
        write_batch_retest_log_file(name, base_dir, primer_df)


def main():
    # Setting up CLI
    parser = ArgumentParser()
    parser.add_argument("--wishlist_file",
                        type=str,
                        help="Path to txt file with gene or guide wishlist. Note: if guides pass batch_retest flag.")
    parser.add_argument("--name",
                        type=str,
                        help="Your name.")
    parser.add_argument("--ai",
                        type=str,
                        default="i",
                        help="i for CRISPRi or a for CRISPRa.")
    parser.add_argument("--guides_per_gene",
                        type=int,
                        default=5,
                        help="Number of guides to order for each gene.")
    parser.add_argument("--order_format",
                        type=str,
                        default='single',
                        help="order types: arrayed, pooled, or single.")
    parser.add_argument("--check_db",
                        action="store_true",
                        help="Whether to check the database for already cloned guides.")
    parser.add_argument("--organism",
                        type=str,
                        default="human",
                        help="Guides that target the mouse or human genome.")
    parser.add_argument("--ntc_frac",
                        type=float,
                        default=0.2,
                        help="The ratio of your library to order as NTCs: 0.2 equates to 20 percent the length of the library."
                        )
    args = parser.parse_args()
    # Reading wishlist genes into list
    file = args.wishlist_file

    if args.order_format.lower() == "pooled":
        gene_list, left_primers, right_primers, lib_num = read_gene_list_pooled(file=file)
        guide_list = []
        primer_df = pd.DataFrame({'gene_symbol': gene_list,
                                  'left_primers': left_primers,
                                  'right_primers': right_primers,
                                  'lib_num': lib_num})
    elif args.order_format.lower() == "batch-retest":
        guide_list, left_primers, right_primers, lib_num = read_gene_list_pooled(file=file)
        gene_list = []
        primer_df = pd.DataFrame({'guide_id': guide_list,
                                  'left_primers': left_primers,
                                  'right_primers': right_primers,
                                  'lib_num': lib_num})
    else:
        guide_list, gene_list = read_wishlist(file=file)
        primer_df = None
    # Defining organism
    organism_name = args.organism
    organism_name = organism_name.lower()
    # Defining base directory to save files to
    if platform == "linux" or platform == "linux2" or platform == "darwin":
        base_dir = "/".join(file.split("/")[:-1])
    elif platform == "win32":
        base_dir = "\\".join(file.split("\\")[:-1])
    else:
        raise Exception("Operating system not recognized.")
    # Writing csv files for ordering guides
    assert args.order_format in ["single", "pooled", "arrayed", "batch-retest"], "Only single, pooled, batch-retest, or arrayed order formats."
    if args.order_format == "pooled" or args.order_format == "batch-retest":
        if args.order_format == "pooled" or args.order_format == "batch-retest":
            assert primer_df is not None, "You must run pyguide-collate and generate primers for pooled ordering."
        else:
            assert primer_df is not None, "You must run pyguide-batch-retest and generate primers for pooled ordering."
        ntc_frac = args.ntc_frac
        order_guides(guide_list,
                     gene_list,
                     name=args.name,
                     ai_status=args.ai,
                     guides_per_gene=args.guides_per_gene,
                     order_format=args.order_format,
                     base_dir=base_dir,
                     check_db=args.check_db,
                     organism=organism_name,
                     primer_df=primer_df,
                     ntc_frac=ntc_frac)
    else:
        order_guides(guide_list,
                     gene_list,
                     name=args.name,
                     ai_status=args.ai,
                     guides_per_gene=args.guides_per_gene,
                     order_format=args.order_format,
                     base_dir=base_dir,
                     check_db=args.check_db,
                     organism=organism_name)


if __name__ == "__main__":
    main()
