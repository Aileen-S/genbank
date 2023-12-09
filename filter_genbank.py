# Taxonomy filter for genbank files

import argparse
from Bio import SeqIO

# Argument parser
# Add option to find only mito genes, or only selected genes.
parser = argparse.ArgumentParser(description="Search GenBank, retrieve gene sequences and save as fasta.")
parser.add_argument('-i', '--input', type=str, help="Input genbank")
parser.add_argument('-t', '--taxonomy', type=str, help="Taxon to extract")
parser.add_argument('-o', '--output', type=str, help='Output genbank')

args = parser.parse_args()         # Process input args from command line

output = open(args.output, 'w')
records = (r for r in SeqIO.parse(args.input, "gb") if args.taxonomy in r.annotations["taxonomy"])
count = SeqIO.write(records, args.output, "fasta")

print(f"{count} records from {args.input} saved to {args.output}")

