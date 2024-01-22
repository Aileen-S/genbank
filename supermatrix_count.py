import argparse
import csv
import re
from Bio import SeqIO

# Argument parser
parser = argparse.ArgumentParser(description="Count number of genes and sequence length in supermatrix")
parser.add_argument("-f", "--fasta", type=str, help="Input supermatrix")
parser.add_argument("-p", "--partitions", type=str, help="Input partition file from catfasta2phyml")

args = parser.parse_args()

# Read partition file
with open(args.partitions) as file:
    parts = {}
    genes = []
    lines = file.readlines()
    for line in lines:
        line = line.strip()
        line = re.split(r'/|\.| = |-', line)
        # Add to dict list for each gene [start pos, end pos, length]
        start = int(line[3]) - 1
        end = int(line[4]) - 1
        parts[line[1]] = [start, end, end - start]
        genes.append(line[1])
print(f'{len(genes)} genes in partition file')

# Write CSV metadata file
with open("count_genes.csv", "w") as file:     # Open output file
    writer = csv.writer(file)               # Name writer object
    writer.writerow(["Taxon", "Count"] + genes)  # Write column names
    records = SeqIO.parse(args.fasta, "fasta")
    x = 0
    for rec in records:
        x += 1
        row = [rec.id, '']
        total = 0
        for gene in genes:
            # Count gaps in each gene
            gaps = rec.seq.count('-', parts[gene][0], parts[gene][1])
            length = parts[gene][2] - gaps
            row.append(length)
            total = total + length
        row[1] = total
        writer.writerow(row)
print(f'{x} taxa written to count_genes.csv')
