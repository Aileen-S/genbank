# Seach GenBank and return number of records for specified taxon and gene(s).
# Only set up for mitochondrial genes at the moment.
# Input taxon, gene and email address.

import urllib
import argparse
from Bio import Entrez

def loadnamevariants():
    conversion = {}
    url = "https://raw.githubusercontent.com/tjcreedy/constants/master/gene_name_variants.txt"
    fullparse = {}
    alltypes = set()
    for line in urllib.request.urlopen(url):
        line = line.decode('utf-8').strip()
        description, variants = line.split(":")
        name, annotype, fullname = description.split(";")
        variants = variants.split(',')
        variants.extend([name, fullname.upper()])

        fullvariants = []
        for v in [name] + variants:
            for g in ['', ' ']:
                v = v.replace(g, '')
                for s in ['', ' GENE', ' '+annotype.upper()]:
                    fullvariants.append(v+s)
                    conversion[v+s] = name

        alltypes.add(annotype)
        fullparse[name] = {'type': annotype, 'variants': fullvariants, 'product': fullname}
    return conversion, alltypes, fullparse

# Argument parser
parser = argparse.ArgumentParser(description="`Count number of GenBank records for specified taxon and gene(s).")
parser.add_argument("-t", "--taxon", type=str, help="Taxon of interest")
parser.add_argument("-g", "--gene", type=str, help="Gene(s) of interest. For multiple genes, input as gene1,gene2,gene3")
parser.add_argument("-e", "--email", type=str, help="Your email address. Required by NCBI.")
args = parser.parse_args()

args = parser.parse_args("-t Agabus -g COX1,ND2,CYTB -e thomas.creedy@gmail.com".split(' '))
Entrez.email = args.email
genes = args.gene.upper()       # Change input genes to upper case
genes = genes.split(",")        # Split genes string into a list


geneconv, types, genefull = loadnamevariants()
stdgenes = [geneconv[g] for g in genes]
names = {sg:genefull[sg]['variants'] for sg in stdgenes}

# Nice work reusing the code from loadnamevariants. Here's a suggestion using an output
# already available from the loadnamevariants function. First get the three outputs,
# then standardise the input list, then just pull out from the full variants dict of dicts
# the list of variants for each gene.


# names = {}                      # Make empty dict for name variants
# for gene in genes:
#     url = "https://raw.githubusercontent.com/tjcreedy/constants/master/gene_name_variants.txt"
#     for line in urllib.request.urlopen(url):
#         line = line.decode('utf-8').strip()
#         description, variants = line.split(":")             # Variants = everything after colon
#         name, annotype, fullname = description.split(";")   # Name = name before first semicolon
#         variants = variants.split(',')                      # Save variants as a list
#         if gene in variants:
#             names[name] = " OR " .join(variants)    # Add to dict, gene as key, variants in string form as value

print(f"Searching {args.taxon} records for genes: {list(names.keys())}")

for gene in names:      # Iterate through dict, new search for each gene names list
    # gene = list(names.keys())[0]
    handle = Entrez.esearch(db="nucleotide", term=f"{args.taxon} AND ({' OR '.join(names[gene])})")
    record = Entrez.read(handle)
    count = int(record["Count"])
    print(str(record["Count"]) + " records found for " + gene)