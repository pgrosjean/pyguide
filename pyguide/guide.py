import pandas as pd
import numpy as np
import os
from typing import List, Dict, Tuple, Union
from argparse import ArgumentParser
from datetime import datetime
from sys import platform
import mygene


####################
# Defining Functions
####################

####################
# I/O Functionality
####################

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
    assert len(df.columns) == 4, "Must run pyguide-collate on gene wishlists before pooled ordering."
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


def filter_database(df: pd.DataFrame, gene_list: List[str]) -> pd.DataFrame:
    """
    This function filters the database.

    Parameters
    ----------
    df : pd.DataFrame
        The dataframe for filtering.
    gene_list : List[str]
        The list of genes used for filtering.

    Returns
    -------
    filtered_df : pd.DataFrame
        The filtered dataframe.
    """
    filtered_df = df.loc[df['gene'].isin(gene_list)]
    return filtered_df


def get_guide_db(gene_list: List[str],
                 ai_status: str,
                 organism: str) -> pd.DataFrame:
    """
    This function reads in the human guide database and filters for genes of interest from the user defined gene list.

    Parameters
    ----------
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
            file = os.path.join(file_path, "CRISPRa_v2.txt")
        else:
            raise Exception("ai status must be either a or i")
    db_df = open_txt_database(file, header=False)  # Database DataFrame
    # Adding column names to database
    db_df = db_df.rename(columns={0: 'name', 1: 'gene', 2: 'score', 3: 'seq'})
    # Filtering full database for genes of interest
    db_df = filter_database(db_df, gene_list)
    return db_df


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
    kl_df = filter_database(kl_df, gene_list)
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
                     base_dir: str):
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
    date = get_current_date()
    filename = get_unique_filename(base_dir, f"order_pooled_{name}_{date}.txt")
    if platform == "linux" or platform == "linux2" or platform == "darwin":
        file = f"{base_dir}/{filename}"
    elif platform == "win32":
        file = f"{base_dir}\\{filename}"
    else:
        file = f"{base_dir}/{filename}"  # Always baseline assumes Linux system
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
            f.write(f"{library_number}\t{left_primer}\t{right_primer}")


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
        print(gene_map['out'])
        return None
    else:
        gene_df = gene_map['out']
        if 'notfound' in gene_df.columns:
            gene_df = gene_df.loc[gene_df['notfound'] != True, :]
        return gene_df


####################
# Main entry point for all guide ordering
####################

def order_guides(gene_list: List[str],
                 name: str,
                 ai_status: str,
                 guides_per_gene: int,
                 order_format: str,
                 base_dir: str,
                 kampmann_lab: bool,
                 organism: str,
                 primer_df: Union[pd.DataFrame, None] = None):
    """
    This function orders guides for cloning with a one-at-a-time method.

    Parameters
    ----------
    gene_list: List[str]
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
    kampmann_lab: bool
        Whether the user is in the Kampmann Lab. Controls access to DataBase.
    organism: str
        Whether to get guides that target the mouse or human.
    primer_df: pd.DataFrame
        Data frame containing information about primers for pooled ordering.
    """
    # Making all str inputs for if else statements lower case
    ai_status = ai_status.lower()
    order_format = order_format.lower()
    assert order_format in ['arrayed', 'single', 'pooled'], "Incorrect order format entered."
    assert guides_per_gene <= 10, "No more than 10 guides per gene allowed."

    # Pulling in sgRNA filtered database
    db_df = get_guide_db(gene_list, ai_status=ai_status, organism=organism)

    # Logging any genes in wishlist not found in database
    missing_genes = [gene for gene in gene_list if gene not in set(db_df['gene'].unique())]
    if len(missing_genes) > 0:
        # Checking missing genes against standardized web API for known genes
        missing_df = query_genes(missing_genes, species=organism)
        if missing_df is not None:
            query_map = make_query_map(missing_df)
            alias_names = np.hstack(list(query_map.keys()))
            symbols = np.hstack(missing_df['symbol'].values)
            missing_filter = list(np.hstack([alias_names, symbols]))
            missing_db_df = get_guide_db(missing_filter, ai_status=ai_status, organism=organism)
            # ^ Returns any of aliases or symbols. So then I need to map these outputs to the original gene name
            # that was requested by the user
            db_df = pd.concat([missing_db_df, db_df])
            mapped_missed_genes = [query_map[gene] for gene in missing_db_df['gene'].values]
            all_mapped_genes = np.hstack([mapped_missed_genes, db_df['gene'].values])
            missing_genes = [gene for gene in gene_list if gene not in all_mapped_genes]
            # Capturing modified queries for returning to user
            mod_query_map = {}  # maps original query to updated name
            for name in missing_db_df['gene']:
                if name in set(list(query_map.keys())):
                    mod_query_map[query_map[name]] = name
        else:
            mod_query_map = None
    else:
        mod_query_map = None
    # Updating gene_list
    gene_list = list(db_df['gene'].values)

    # Collate guides
    collated_df = collate_guides(guides_per_gene, db_df)

    # Single Guide Ordering
    if order_format == "single":
        if kampmann_lab is True and organism == "human":
            # Pulling in local Kampmann Lab Database
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
        write_basic_log_file(missing_genes,
                             mod_query_map,
                             name,
                             base_dir,
                             "arrayed")

    # Pooled guide ordering
    elif order_format == "pooled":
        assert primer_df is not None, "must provide primer information when pooled ordering."
        # changing gene names in the primer df if necessary
        if mod_query_map is not None:
            primer_df['gene_symbol'] = primer_df['gene_symbol'].apply(lambda x: mod_query_map.get(x, x))
        collated_df = collated_df.merge(primer_df, left_on='gene', right_on='gene_symbol')
        # Writing the pooled ordering txt file
        write_pooled_txt(collated_df, name, base_dir)
        write_pooled_log_file(missing_genes, mod_query_map, name, base_dir, primer_df)


def main():
    ## TODO: Get rid of the Kampmann Lab flag and replace with the AWS RDS Queries
    # Setting up CLI
    parser = ArgumentParser()
    parser.add_argument("--wishlist_file",
                        type=str,
                        help="Path to txt file with gene wishlist.")
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
    parser.add_argument("--kampmann_lab",
                        action="store_true",
                        help="Are you a member of the Kampmann Lab?")
    parser.add_argument("--organism",
                        type=str,
                        default="human",
                        help="Guides that target the mouse or human genome.")
    args = parser.parse_args()
    # Reading wishlist genes into list
    file = args.wishlist_file

    if args.order_format.lower() == "pooled":
        gene_list, left_primers, right_primers, lib_num = read_gene_list_pooled(file=file)
        primer_df = pd.DataFrame({'gene_symbol': gene_list,
                                  'left_primers': left_primers,
                                  'right_primers': right_primers,
                                  'lib_num': lib_num})
    else:
        gene_list = read_gene_list(file=file)
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
    assert args.order_format in ["single", "pooled", "arrayed"], "Only single, pooled, or arrayed order formats."
    if args.order_format == "pooled":
        assert primer_df is not None, "You must run pyguide-collate and generate primers for pooled ordering."
        order_guides(gene_list,
                     name=args.name,
                     ai_status=args.ai,
                     guides_per_gene=args.guides_per_gene,
                     order_format=args.order_format,
                     base_dir=base_dir,
                     kampmann_lab=args.kampmann_lab,
                     organism=organism_name,
                     primer_df=primer_df)
    else:
        order_guides(gene_list,
                     name=args.name,
                     ai_status=args.ai,
                     guides_per_gene=args.guides_per_gene,
                     order_format=args.order_format,
                     base_dir=base_dir,
                     kampmann_lab=args.kampmann_lab,
                     organism=organism_name)


if __name__ == "__main__":
    main()
