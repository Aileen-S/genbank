
#python3 get_concat_recs.py -e mixedupvoyage@gmail.com -t Megadytes

import argparse
from Bio import SeqIO


# Argument parser
# Add option to find only mito genes, or only selected genes.
parser = argparse.ArgumentParser(description="Remove duplicate sequences from fasta.")
parser.add_argument("-i", "--input", type=str, help="Input fasta")
parser.add_argument('-o', '--output', type=str, help="Output fasta with duplicates removed")
parser.add_argument('-d', '--dups', type=str, help="File for removed duplicate sequences")

args = parser.parse_args()         # Process input args from command line
#args = argparse.Namespace(input='ND1.fasta')

sequences = []
check = set()
seq_ids = []
dup_ids = []
with open(args.input) as file:
    record = SeqIO.parse(file, "fasta")
    for rec in record:
        if rec.seq not in sequences:
            sequences.append(rec.seq)
            check.add(rec.seq)
            seq_ids.append(rec.id)
        else:
            dup_ids.append(rec.id)

records = (r for r in SeqIO.parse(args.input, "fasta") if r.id in seq_ids)
count = SeqIO.write(records, args.output, "fasta")
print(f"Saved {count} records from {args.input} to {args.output}")

if args.dups:
    dups = (r for r in SeqIO.parse(args.input, "fasta") if r.id in dup_ids)
    count = SeqIO.write(dups, args.dups, "fasta")
    print(f"{count} duplicate records removed and saved to {args.dups}")




# Set record length as 0, iterate through records and replace whenever another sequence is longer.
