![example workflow](https://github.com/pgrosjean/pyguide/actions/workflows/python_package-linux.yml/badge.svg)
![example workflow](https://github.com/pgrosjean/pyguide/actions/workflows/python_package-macos.yml/badge.svg)
![example workflow](https://github.com/pgrosjean/pyguide/actions/workflows/python_package-windows.yml/badge.svg)

# PYGUIDE
Tools for ordering gRNA and maintaining gRNA libraries for CRISPRi/a work.

# Installation
```bash
git clone github.com/pgrosjean/pyguide
cd pyguide
pip install -e .
```

# Usage
## Commandline Usage
```bash
# view options for ordering guides
pyguide-order --help
```

```bash
# EXAMPLE: Ordering 5 CRISPRi human guides per gene for one at a time cloning
pyguide-order --wishlist_file=/path/to/wishlist_file --name=Name --ai=i --guides_per_gene=5 --order_format=single

# EXAMPLE: Ordering 3 CRISPRa human guides per gene for arrayed library
pyguide-order --wishlist_file=/path/to/wishlist_file --name=Name --ai=a --guides_per_gene=3 --order_format=arrayed

# EXAMPLE: Ordering 2 CRISPRi human guides per gene for pooled library
pyguide-order --wishlist_file=/path/to/wishlist_file --name=Name --ai=i --guides_per_gene=2 --order_format=pooled

# EXAMPLE: Ordering 3 CRISPRi mouse guides per gene for one at a time cloning
pyguide-order --wishlist_file=/path/to/wishlist_file --name=Name --ai=i --guides_per_gene=3 --order_format=single --organism=mouse
```

### Flags for gRNA ordering using pyguide-order

--wishlist_file : This is the file that contains a list of your target genes (see testing/examples/gene_list.txt for example)

--name : Your name

--ai : Whether to use CRISPRa v2 vs CRISPRi v2 library (pass 'a' or 'i' for this flag)

--guides_per_gene : Number of guides to order per gene

--order_gormat : Whether to order single guides (IDT), arrayed library (IDT), or pooled library (Agilent) (pass 'single', 'arrayed', or 'pooled' for this flag)

--organism: Whether to generate guides for the human or mouse genome (pass 'mouse' or 'human' to this flag) defaults to human.

--kampmann_lab: Whether you are part of the Kampmann lab and have the database configured. Only pass true if you are in the lab.

### pyguide-order outputs
Upon running pyguide-order two files will be saved to the same directory as the wishlist file. 

#### (1) Log File

- **For 'arrayed' and 'pooled' order formats:**
This file tells you if any of the genes that were on your wishlist could not be found in the guide databse. 

- **For 'single' order format:**
In addition to telling you if a requested gene is not in the databse, this log file will also tell you where to find guides that have already been cloned. Note the guides that have already been cloned will not be in your ordering file.


#### (2) Order File(s)

- **For 'single' order format:**
This file will be a csv file that starts with order_single_ and can be uploaded directly to the [IDT OligoEntry](https://www.idtdna.com/site/order/oligoentry) Bulk Input for ordering guides for cloning in a one at a time manner.

- **For 'arrayed' order format:**
There will be at least two files produced for this order format, all of which will be csv files. There will always be at least two plate corresponding to the top and bottom oligo pairs, which will match to make for easy arrayed guide library prep. If there are are more than 96 guides that you are ordering this will expand to more than one pair of plates. These files will begin with order_arrayed_ and end with plate_#.csv. These files can be uploaded to the [IDT DNA Plates order page](https://www.idtdna.com/site/order/plate/index/dna/1799)

- **For 'Pooled' order format:**
This file will be a text file that starts with order_pooled_ and can be used to order a pooled library from Agilent.
