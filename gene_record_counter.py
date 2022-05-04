# Seach GenBank and return number of records for specified taxon and gene(s).
# Input taxon, gene and email address.

import urllib
import argparse
from Bio import Entrez

# Argument parser
parser = argparse.ArgumentParser(description="`Count number of GenBank records for specified taxon and gene(s).")
parser.add_argument("-t", "--taxon", type=str, help="Taxon of interest")
parser.add_argument("-g", "--gene", type=str, help="Gene(s) of interest. For multiple genes, input as gene1,gene2,gene3")
parser.add_argument("-e", "--email", type=str, help="Your email address. Required by NCBI.")
args = parser.parse_args()

Entrez.email = args.email
genes = args.gene.split(",")    # Split genes string into a list
names = {}                      # Make empty dict for name variants
for gene in genes:
    url = "https://raw.githubusercontent.com/tjcreedy/constants/master/gene_name_variants.txt"
    for line in urllib.request.urlopen(url):
        line = line.decode('utf-8').strip()
        description, variants = line.split(":")
        name, annotype, fullname = description.split(";")
        variants = variants.split(',')
        variants.extend([name, fullname.upper()])
        if gene in variants:
            names[gene] = " OR " .join(variants)    # Add to dict, gene as key, variants in string form as value

for gene in names:      # Iterate through dict, new search for each gene names list
    handle = Entrez.esearch(db="nucleotide", term=f"{args.taxon} AND ({names[gene]})")
    record = Entrez.read(handle)
    count = int(record["Count"])
    print(str(record["Count"]) + " records found for " + gene)

