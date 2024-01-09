from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Try each reading frame for sequences which did not align")
parser.add_argument("-i", "--input", type=str, help="Input fasta to be filtered")
parser.add_argument("-f", "--found", type=str, help="Output fasta for ids in search file")
parser.add_argument("-n", "--notfound", type=str, help="Output fasta for ids not in search file")
parser.add_argument("-s", "--search", type=str, help="File with IDs to search for")

args = parser.parse_args()         # Process input args from command line

search = []
with open(args.search) as infile:
    lines = infile.readlines()
    for line in lines:
        search.append(line.strip())

found = {}
notfound = {}
with open(args.input) as file:
    records = (r for r in SeqIO.parse(args.input, "fasta") if r.id in search)
    for rec in records:
        found[rec.id[:-1]] = rec.seq
    records = (r for r in SeqIO.parse(args.input, "fasta") if r.id not in search)
    for rec in records:
        notfound[rec.id[:-1]] = rec.seq

if args.found:
    with open(args.found, "w") as outfile:
        for rec_id, seq in found.items():
            outfile.write(f'{rec_id}1\n{seq}\n')
            outfile.write(f'{rec_id}2\n{seq}\n')
            outfile.write(f'{rec_id}3\n{seq}\n')


if args.notfound:
    with open(args.notfound, "w") as outfile:
        for rec_id, seq in notfound.items():
            outfile.write(f'{rec_id}1\n{seq}\n')
            outfile.write(f'{rec_id}2\n{seq}\n')
            outfile.write(f'{rec_id}3\n{seq}\n')
