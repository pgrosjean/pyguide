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
### Flags for gRNA ordering using pyguide-order

--wishlist_file : This is the file that contains a list of your target genes (see testing/examples/gene_list.txt for example)

--name : Your name

--ai : Whether to use CRISPRa v2 vs CRISPRi v2 library (pass 'a' or 'i' for this flag)

--guides_per_gene : Number of guides to order per gene

--order_gormat : Whether to order single guides (IDT), arrayed library (IDT), or pooled library (Agilent) (pass 'single', 'arrayed', or 'pooled' for this flag)
```bash
# EXAMPLE: Ordering 5 CRISPRi guides per gene for one at a time cloning
pyguide-order --wishlist_file=/path/to/wishlist_file --name=Name --ai=i --guides_per_gene=5 --order_format=single

# EXAMPLE: Ordering 3 CRISPRa guides per gene for arrayed library
pyguide-order --wishlist_file=/path/to/wishlist_file --name=Name --ai=a --guides_per_gene=3 --order_format=arrayed

# EXAMPLE: Ordering 2 CRISPRi guides per gene for pooled library
pyguide-order --wishlist_file=/path/to/wishlist_file --name=Name --ai=i --guides_per_gene=2 --order_format=pooled
```
