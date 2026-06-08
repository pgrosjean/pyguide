![example workflow](https://github.com/pgrosjean/pyguide/actions/workflows/python_package-linux.yml/badge.svg)
![example workflow](https://github.com/pgrosjean/pyguide/actions/workflows/python_package-macos.yml/badge.svg)
![example workflow](https://github.com/pgrosjean/pyguide/actions/workflows/python_package-windows.yml/badge.svg)

<img width="349" alt="image" src="https://github.com/user-attachments/assets/bb598808-3355-4c23-b85a-984abd76701a" />

# PYGUIDE
**Overview**: Tools for ordering gRNA and maintaining gRNA libraries for CRISPRi/a work. All functionalities assume that you are ordering guides for use in the pMK1334 plasmid (which you can request from the Kampmann Lab) and that you are using gRNAs from the libraries defined in the [Horlbeck et al 2016 paper](https://elifesciences.org/articles/19760#content).

# Installation

**Requirements:** Python 3.12+, [uv](https://docs.astral.sh/uv/)

```bash
git clone https://github.com/pgrosjean/pyguide.git
cd pyguide
```

**CLI only** (lighter install — all four `pyguide-*` commands):
```bash
uv sync
```

**Streamlit app** (includes the web UI):
```bash
uv sync --extra app
```

**Development** (CLI + tests):
```bash
uv sync --extra dev
```

**Everything** (CLI + app + tests):
```bash
uv sync --extra app --extra dev
```

**Tiling library + BigWig output** (includes pyBigWig for genome browser tracks):
```bash
uv sync --extra tiling
```

# Usage (Streamlit App)

## Using the Streamlit App

```bash
uv run streamlit run pyguide/app.py
```

The app supports all order formats available from the CLI. Wishlist files are `.txt` files containing either Gene Symbols (e.g. `APOE`) or sgRNA names (e.g. `BIRC7_+_61867189.23-P1P2`) from the libraries in `/data/`. For the **pooled-seq** format, upload tab-delimited sequence files (`name<TAB>20nt spacer`) instead.


# Usage (CLI)
## Using the CLI
```bash
# view options for pyguide-order
pyguide-order --help

# view options for pyguide-collate
pyguide-collate --help

# view options for pyguide-batch-retest
pyguide-batch-retest --help

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
# EXAMPLE: Ordering 5 CRISPRi human guides per gene for single guide format (checking for guides that have already been cloned in /data/human_sgrnas.txt)
pyguide-order --wishlist_file /path/to/wishlist_file --name Name --ai i --guides_per_gene 5 --order_format single --check_db
```

(2) Ordering human CRISPRi/a guides **without** checking for guides that have already been cloned
```bash
# EXAMPLE: Ordering 4 CRISPRa human guides per gene for single guide format
pyguide-order --wishlist_file /path/to/wishlist_file --name Name --ai a --guides_per_gene 4 --order_format single
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
pyguide-order --wishlist_file /path/to/wishlist_file --name Name --ai a --guides_per_gene 3 --order_format arrayed

# EXAMPLE: Ordering 3 CRISPRa mouse guides per gene for arrayed format
pyguide-order --wishlist_file /path/to/wishlist_file --name Name --ai a --guides_per_gene 3 --order_format arrayed --organism mouse
```

### pyguide-order outputs for arrayed ordering
Upon running pyguide-order two files will be saved to the same directory as the wishlist file. 

**(1) Log File**
This file contains information on:
- Who generated the order file and which date it was generated.
- If any genes in your original wishlist file were converted to an alias for guide ordering.
- If any genes in your original wishlist file are not valid genes for ordering guides.

**(2) Order Files**
There will be at least two files produced for this order format, all of which will be csv files. There will always be at least two plates corresponding to the top and bottom oligo pairs, which will match to make for easy arrayed guide cloning. If there are more than 96 guides that you are ordering this will expand to more than one pair of plates. These files will begin with order_arrayed_ and end with plate_#.csv. These files can be uploaded to the [IDT DNA Plates order page](https://www.idtdna.com/site/order/plate/index/dna/1799) for ordering.



















## Pooled Guide Ordering Usage

**Note: For pooled guide ordering you must run pyguide-collate before running pyguide-order.**
This is to enable multiple pooled libraries to be ordered at once that can then be selectively amplified out of the oligo pool using libary specific primers.

### Flags for pyguide-collate

**Flags you must provide**:
- **wishlist_files**: Paths to one or more wishlist files (one wishlist file per pooled library)

**Optional Flags**:
- **primer_file**: Path to a user-defined primer-file, otherwise random pairs of primers will be used from the /data/pooled_primers.txt file. (If you want to specify your own primers then generate a text file following the format of /data/pooled_primers.txt)

### Output of pyguide-collate
A collated pooled wishlist file that contains the genes of interest and the primers that will be used for sequencing will be saved to the same directory as your user-defined wishlist that you provide to pyguide-collate.

### Flags for pyguide-order for pooled ordering

**Flags you must provide**:
- **wishlist_file**: Path to your file containing a list of the genes for which you wish to order guides. (see testing/examples/gene_list.txt for example)
- **name**: Your name.
- **ai**: Whether you want to order for CRISPR interference (i) or activation (a).
- **guides_per_gene**: The number of guides per gene you wish to order.
- **order_format**: Whether to order single guides (IDT), arrayed library (IDT), or pooled library (Agilent). Pass arrayed for arrayed ordering and cloning protocol.

**Optional Flags**:
- **ntc_frac**: The fraction of the library to generate NTCs for (e.g. 0.1 will make enough non-targeting controls to be equal to 10% of the library size)
- **organism**: Whether to order guides that target mouse genes or human genes (defaults to human).


**Example use cases**:
(1) Ordering a single library in a single DNA Oligo Pool with randomly assigned primers.
```bash
# EXAMPLE: Ordering 3 CRISPRi mouse guides per gene for a single pooled library using random primers
## Step 1: Generate collated wishlist file
pyguide-collate --wishlist_files /path/to/wishlist_file
## Step 2: Using collated wishlist file to order the DNA oligo pool
pyguide-order --wishlist_file /path/to/collated_wishlist_file --name Name --ai i --guides_per_gene 3 --order_format pooled --organism mouse
```

(2) Ordering multiple libraries in a single DNA Oligo Pool with randomly assigned primers.
```bash
# EXAMPLE: Ordering 3 CRISPRa human guides per gene for a single pooled library using random primers
## Step 1: Generate collated wishlist file
pyguide-collate --wishlist_files /path/to/wishlist_file_1 /path/to/wishlist_file_2
## Step 2: Using collated wishlist file to order the DNA oligo pool
pyguide-order --wishlist_file /path/to/collated_wishlist_file --name Name --ai a --guides_per_gene 3 --order_format pooled
```

(3) Ordering a single libary in a single DNA Oligo Pool with user-defined primers.
```bash
# EXAMPLE: Ordering 4 CRISPRi human guides per gene for a single pooled library using random primers
## Step 1: Generate collated wishlist file
pyguide-collate --wishlist_files /path/to/wishlist_file --primer_file /path/to/primer/file
## Step 2: Using collated wishlist file to order the DNA oligo pool
pyguide-order --wishlist_file /path/to/collated_wishlist_file --name Name --ai i --guides_per_gene 4 --order_format pooled
```

(4) Ordering multiple libraries in a single DNA Oligo Pool with user-defined primers
```bash
# EXAMPLE: Ordering 3 CRISPRa human guides per gene for single pooled library using random primers
## Step 1: Generate collated wishlist file
pyguide-collate --wishlist_files /path/to/wishlist_file_1 /path/to/wishlist_file_2 --primer_file /path/to/primer/file
## Step 2: Using collated wishlist file to generate a pooled library for 3 CRISPRa human guides per gene.
pyguide-order --wishlist_file /path/to/collated_wishlist_file --name Name --ai a --guides_per_gene 3 --order_format pooled
```

### pyguide-order outputs for pooled ordering
Upon running pyguide-order two files will be saved to the same directory as the wishlist file. 

**(1) Log File**
This file contains information on:
- Who generated the order file and on which date it was generated.
- Which primers are assigned to which libraries.
- If any genes in your original wishlist file were converted to an alias for guide ordering.
- If any genes in your original wishlist file are not valid genes for ordering guides.

**(2) Order Files**
This is a text file that starts with order_pooled_ and can be used to order a DNA oligo pool from Agilent.















## Batch-retest Guide Ordering Usage
**Note: For batch retest guide ordering you must run pyguide-batch-retest before running pyguide-order.**
This is to enable the usage of batch-retest specific primers for ordering a pooled library for batch retest.

### Flags for pyguide-batch-retest
- **wishlist_file**: Paths to a single wishlist file containing sgRNA names (Note: you need to provide specific guide ids not gene names)

**Optional Flags**:
- **primer_file**: Path to a user-defined primer-file, otherwise a random pair of primers will be used from the /data/pooled_primers.txt file. (If you want to specify your own primers then generate a text file following the format of /data/pooled_primers.txt)

### Output of pyguide-batch-retest
A collated pooled wishlist file that contains the genes of interest and the primers that will be used for sequencing will be saved to the same directory as your user-defined wishlist that you provide to pyguide-collate.

### Flags for pyguide-order for batch-retest ordering

**Flags you must provide**:
- **wishlist_file**: Path to your file containing a list of the genes for which you wish to order guides. (see testing/examples/gene_list.txt for example)
- **name**: Your name.
- **ai**: Whether you want to order for CRISPR interference (i) or activation (a).
- **guides_per_gene**: The number of guides per gene you wish to order.
- **order_format**: Whether to order single guides (IDT), arrayed library (IDT), pooled library (Agilent), or batch-retest library (Agilent). Pass arrayed for arrayed ordering and cloning protocol.

**Optional Flags**:
- **ntc_frac**: The fraction of the library to generate NTCs for (e.g. 0.1 will make enough non-targeting controls to be equal to 10% of the library size)
- **organism**: Whether to order guides that target mouse genes or human genes (defaults to human).

### pyguide-order outputs for batch-retest ordering
Upon running pyguide-order two files will be saved to the same directory as the wishlist file. 

**(1) Log File**
This file contains information on:
- Who generated the order file and on which date it was generated.
- Which primers are assigned to which libraries.
- If any genes in your original wishlist file were converted to an alias for guide ordering.
- If any genes in your original wishlist file are not valid genes for ordering guides.

**(2) Order Files**
This is a text file that starts with order_pooled_ and can be used to order a DNA oligo pool from Agilent.

**Example use case**:
(1) Ordering a batch-retest library
```bash
# EXAMPLE: Ordering CRISPRi human guides from a wishlist of guide_ids
## Step 1: Generate batch-retest wishlist file
pyguide-batch-retest --wishlist_file /path/to/batch_retest_wishlist_file_guides
## Step 2: Using batch-retest returned wishlist file to order the DNA oligo pool
pyguide-order --wishlist_file /path/to/batch_retest_generated_wishlist_file --name Name --ai i --order_format batch-retest --organism human
```








## pyguide-check-seq usage

### Flags for pyguide-check-seq

**Flags you must provide**:
- **file_dir**: Path to the directory containing .seq files generated from sanger sequencing of the guides you ordered (see testing/example/seq_files as an example directory).
- **ai**: Whether the guides that you are checking the sequencing results of are CRISPR interference (i) or activation (a).
- **organism**: Whether the guides that you are checking the sequencing results of target mouse genes or human genes.
- **name**: Your username
- **backbone**: The backbone you are using. Either pMK1334 or pLG15.
- **update_db**: Only include this flag if you are checking sequencing for single guides that should be put into the local database.

### pyguide-check-seq outputs
Upon running pyguide-check-seq a text file will be generated with two tab seperated columns, the first column with the .seq file name and the second column with the corresponding guide. If a .seq file does not correspond to a guide then it will not be returned in the text file.

## Genomic Tiling Library Usage

### Prerequisites

**GuideScan2** must be installed separately (not available on PyPI):
```bash
conda install -c bioconda guidescan2
```

**hg38 GuideScan2 index** (~2.2 GB) — download with the provided script:
```bash
bash scripts/download_guidescan_data.sh
# Default saves to ~/.pyguide_data/
# Use --data-dir to choose a custom location
```

**BigWig output** (optional) — install pyBigWig:
```bash
uv sync --extra tiling
```

### Flags for pyguide-tiling

**Required (one of):**
- **--coordinates**: Single coordinate region, e.g. `chr1:12345-12395` (1-based, closed)
- **--coordinates_file**: Path to a file with one coordinate region per line

**Required:**
- **--index**: Path to the GuideScan2 hg38 index file downloaded above
- **--output**: Path for the output tab-delimited sequence file

**Optional:**
- **--specificity**: Minimum GuideScan2 specificity score (default: 0.2)
- **--hamming**: Minimum pairwise Hamming distance between any two selected guides (default: 4)
- **--guides_per_region**: Maximum number of guides per region, evenly spread across positions
- **--bigwig**: Also write a BigWig coverage track (`.bw` file). Requires `uv sync --extra tiling`.

### Example use cases

**(1) Tile a single promoter region:**
```bash
pyguide-tiling \
  --index ~/.pyguide_data/hg38.index \
  --coordinates chr1:12345-12595 \
  --output tiling_guides.txt
```

**(2) Tile multiple regions from a file:**
```bash
# regions.txt contains one region per line:
# chr1:12345-12595
# chr7:117548601-117548801

pyguide-tiling \
  --index ~/.pyguide_data/hg38.index \
  --coordinates_file regions.txt \
  --output tiling_guides.txt \
  --guides_per_region 20
```

**(3) Tile with BigWig coverage output:**
```bash
pyguide-tiling \
  --index ~/.pyguide_data/hg38.index \
  --coordinates chr7:117548601-117548801 \
  --output tiling_guides.txt \
  --bigwig
# Produces: tiling_guides.txt and tiling_guides.bw
```

### pyguide-tiling output

A tab-delimited file with two columns: `name<TAB>20nt_spacer`. Guide names encode genomic position and strand (e.g. `chr1:12350_fwd`, `chr1:12370_rev`).

This file feeds directly into the pooled-seq ordering pipeline:
```bash
# Step 1: generate tiling guides
pyguide-tiling --index ~/.pyguide_data/hg38.index \
  --coordinates chr1:12345-12595 --output tiling_guides.txt

# Step 2: collate (assigns library primers)
pyguide-collate-seq --sequence_files tiling_guides.txt

# Step 3: order the pooled oligo library
pyguide-order --wishlist_file collated_seq_wishlist_YY_MM_DD.txt \
  --name Name --ai i --guides_per_gene 1 --order_format pooled-seq
```

