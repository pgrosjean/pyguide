![example workflow](https://github.com/pgrosjean/pyguide/actions/workflows/python_package-linux.yml/badge.svg)
![example workflow](https://github.com/pgrosjean/pyguide/actions/workflows/python_package-macos.yml/badge.svg)
![example workflow](https://github.com/pgrosjean/pyguide/actions/workflows/python_package-windows.yml/badge.svg)

# PYGUIDE
**Overview**: Tools for ordering gRNA and maintaining gRNA libraries for CRISPRi/a work. All functionalities assume that you are ordering guides for use in the pMK1334 plasmid (which you can request from the Kampmann Lab) and that you are using gRNAs from the libraries defined in the [Horlbeck et al 2016 paper](https://elifesciences.org/articles/19760#content).

**Database functionality**: The current implementation for checking if gRNAs have already been cloned in this current repository uses a user-provided CSV file /data/human_sgrnas.csv. **Note: this is only implemented for single order format for the human species.** Current work is being done to implement AWS databse functionality. Stay tuned...

# Installation
```bash
git clone github.com/pgrosjean/pyguide
cd pyguide
pip install -e .
```









# Usage
## Overall Commandline Usage
```bash
# view options for pyguide-order
pyguide-order --help

# view options for pyguide-collate
pyguide-collate --help

# view options for pyguide-check-seq
pyguide-check-seq --help
```











## Single Guide Ordering Usage

### Flags for pyguide-order

**Flags you must provide**:
- **wishlist_file**: Path to your file containing a list of the genes for which you wish to order guides. (see testing/examples/gene_list.txt for example)
- **name**: Your name.
- **ai**: Whether you want to order for CRISPR interference (i) or activation (a).
- **guides_per_gene**: The number of guides per gene you wish to order.
- **order_format**: Whether to order single guides (IDT), arrayed library (IDT), or pooled library (Agilent). Pass single for single guide ordering and cloning protocol.

**Optional Flags**:
- **organism**: Whether to order guides that target mouse genes or human genes (defaults to human).
- **check_db**: Pass this flag if you wish to check a databse for guides that have already been cloned.

**Example use cases**:

(1) Ordering human CRISPRi/a guides and checking database for guides that have already been cloned
```bash
# EXAMPLE: Ordering 5 CRISPRi human guides per gene for single guide format (checking for guides that have already been cloned)
pyguide-order --wishlist_file=/path/to/wishlist_file --name=Name --ai=i --guides_per_gene=5 --order_format=single --check_db
```

(2) Ordering human CRISPRi/a guides **without** checking for guides that have already been cloned
```bash
# EXAMPLE: Ordering 4 CRISPRa human guides per gene for single guide format
pyguide-order --wishlist_file=/path/to/wishlist_file --name=Name --ai=a --guides_per_gene=4 --order_format=single
```

### pyguide-order outputs for single guide ordering
Upon running pyguide-order two files will be saved to the same directory as the wishlist file. 

**(1) Log File**
This file contains information on:
- Who generated the order file and which date it was generated.
- If any genes in your original wishlist file were converted to an alias for guide ordering.
- If any genes in your original wishlist file are not valid genes for ordering guides.
- If you passed the check_db flag then this will also tell you if any of the guides you requsted have already been ordered and if so what their corresponding guide number is.

**(2) Order File**
This file will be a csv file that starts with order_single_ and can be uploaded directly to the [IDT OligoEntry](https://www.idtdna.com/site/order/oligoentry) Bulk Input for ordering guides for cloning in a one at a time manner.















## Arrayed Guide Ordering Usage

### Flags for pyguide-order

**Flags you must provide**:
- **wishlist_file**: Path to your file containing a list of the genes for which you wish to order guides. (see testing/examples/gene_list.txt for example)
- **name**: Your name.
- **ai**: Whether you want to order for CRISPR interference (i) or activation (a).
- **guides_per_gene**: The number of guides per gene you wish to order.
- **order_format**: Whether to order single guides (IDT), arrayed library (IDT), or pooled library (Agilent). Pass arrayed for arrayed ordering and cloning protocol.

**Optional Flags**:
- **organism**: Whether to order guides that target mouse genes or human genes (defaults to human).

**Example use case**:
(1) Ordering an arrayed library of any size
```bash
# EXAMPLE: Ordering 3 CRISPRa human guides per gene for arrayed format
pyguide-order --wishlist_file=/path/to/wishlist_file --name=Name --ai=a --guides_per_gene=3 --order_format=arrayed

# EXAMPLE: Ordering 3 CRISPRa mouse guides per gene for arrayed format
pyguide-order --wishlist_file=/path/to/wishlist_file --name=Name --ai=a --guides_per_gene=3 --order_format=arrayed --organism=mouse
```

### pyguide-order outputs for arrayed ordering
Upon running pyguide-order two files will be saved to the same directory as the wishlist file. 

**(1) Log File**
This file contains information on:
- Who generated the order file and which date it was generated.
- If any genes in your original wishlist file were converted to an alias for guide ordering.
- If any genes in your original wishlist file are not valid genes for ordering guides.

**(2) Order Files**
There will be at least two files produced for this order format, all of which will be csv files. There will always be at least two plates corresponding to the top and bottom oligo pairs, which will match to make for easy arrayed guide cloning. If there are are more than 96 guides that you are ordering this will expand to more than one pair of plates. These files will begin with order_arrayed_ and end with plate_#.csv. These files can be uploaded to the [IDT DNA Plates order page](https://www.idtdna.com/site/order/plate/index/dna/1799) for ordering.



















## Pooled Guide Ordering Usage

**Note: For pooled guide ordering you must run pyguide-collate before running pyguide-order.**
This is to enable multiple pooled libraries to be ordered at once that can then be selectively amplified out of the oligo pool using libary specific primers.

### Flags for pyguide-collate

**Flags you must provide**:
- **wishlist_files**: Paths to one or more wishlist files (one wishlist file per pooled library)

**Optional Flags**:
- **primer_file**: Path to a user-defined primer-file, otherwise random pairs of primers will be used from the /data/pooled_primers.txt file. (If you want to specify your own primers then generate a text file following the format of /data/pooled_primers.txt)

### Output of pyguide-collate
A collated pooled wishlist file that contains the genes of interested and the primers that will be used for sequencing will be saved to the same directory as your user-defined wishlist that you provide to pyguide-collate.

### Flags for pyguide-order for pooled ordering

**Flags you must provide**:
- **wishlist_file**: Path to your file containing a list of the genes for which you wish to order guides. (see testing/examples/gene_list.txt for example)
- **name**: Your name.
- **ai**: Whether you want to order for CRISPR interference (i) or activation (a).
- **guides_per_gene**: The number of guides per gene you wish to order.
- **order_format**: Whether to order single guides (IDT), arrayed library (IDT), or pooled library (Agilent). Pass arrayed for arrayed ordering and cloning protocol.

**Optional Flags**:
- **organism**: Whether to order guides that target mouse genes or human genes (defaults to human).


**Example use cases**:
(1) Ordering a single libary in a single DNA Oligo Pool with randomly assigned primers.
```bash
# EXAMPLE: Ordering 3 CRISPRi mouse guides per gene for a single pooled library using random primers
## Step 1: Generate collated wishlist file
pyguide-collate --wishlist_files=/path/to/wishlist_file
## Step 2: Using collated wishlist file to order the DNA oligo pool
pyguide-order --wishlist_file=/path/to/collated_wishlist_file --name=Name --ai=i --guides_per_gene=3 --order_format=pooled --organism=mouse
```

(2) Ordering multiple libraries in a single DNA Oligo Pool with randomly assigned primers.
```bash
# EXAMPLE: Ordering 3 CRISPRa human guides per gene for single pooled library using random primers
## Step 1: Generate collated wishlist file
pyguide-collate --wishlist_files /path/to/wishlist_file_1 /path/to/wishlist_file_2
## Step 2: Using collated wishlist file to order the DNA oligo pool
pyguide-order --wishlist_file=/path/to/collated_wishlist_file --name=Name --ai=a --guides_per_gene=3 --order_format=pooled
```

(3) Ordering a single libary in a single DNA Oligo Pool with user-defined primers.
```bash
# EXAMPLE: Ordering 4 CRISPRi human guides per gene for a single pooled library using random primers
## Step 1: Generate collated wishlist file
pyguide-collate --wishlist_files=/path/to/wishlist_file --primer_file=/path/to/primer/file
## Step 2: Using collated wishlist file to order the DNA oligo pool
pyguide-order --wishlist_file=/path/to/collated_wishlist_file --name=Name --ai=i --guides_per_gene=4 --order_format=pooled
```

(4) Ordering multiple libraries in a single DNA Oligo Pool with user-defined primers
```bash
# EXAMPLE: Ordering 3 CRISPRa human guides per gene for single pooled library using random primers
## Step 1: Generate collated wishlist file
pyguide-collate --wishlist_files /path/to/wishlist_file_1 /path/to/wishlist_file_2 --primer_file /path/to/primer/file
## Step 2: Using collated wishlist file to generate a pooled library for 3 CRISPRa human guides per gene.
pyguide-order --wishlist_file=/path/to/collated_wishlist_file --name=Name --ai=a --guides_per_gene=3 --order_format=pooled
```

### pyguide-order outputs for pooled ordering
Upon running pyguide-order two files will be saved to the same directory as the wishlist file. 

**(1) Log File**
This file contains information on:
- Who generated the order file and which date it was generated.
- Which primers are assigned to which libraries.
- If any genes in your original wishlist file were converted to an alias for guide ordering.
- If any genes in your original wishlist file are not valid genes for ordering guides.

**(2) Order Files**
This is a text file that starts with order_pooled_ and can be used to order a DNA oligo pool from Agilent.











## pyguide-check-seq usage

### Flags for pyguide-check-seq

**Flags you must provide**:
- **file_dir**: Path to the directory containing .seq files generated from sanger sequencing of the guides you ordered (see testing/example/seq_files as an example directory).
- **ai**: Whether the guides that you are checking the sequencing results of are CRISPR interference (i) or activation (a).
- **organism**: Whether the guides that you are checking the sequencing results of target mouse genes or human genes.

### pyguide-check-seq outputs
Upon running pyguide-check-seq a text file will be generated with two tab seperated columns, the first column with the .seq file name and the second column with the corresponding guide. If a .seq file does not correspond to a guide then it will not be returned in the text file.

