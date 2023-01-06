from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Filter fasta file from list of IDs.")
parser.add_argument("-i", "--input", type=str, help="File to be filtered")
parser.add_argument("-f", "--filter", type=str, help="ID list")
parser.add_argument("-o", "--output", type=str, help="Output file")
args = parser.parse_args()

with open(args.filter) as id_handle:
    wanted = set(line.rstrip("\n").split(None, 1)[0] for line in id_handle)

print(f"Found {len(wanted)} unique identifiers in {args.filter}.")
print(f"Searching {args.input}.")

records = (r for r in SeqIO.parse(args.input, "fasta") if r.id in wanted)
count = SeqIO.write(records, args.output, "fasta")

if count < len(wanted):
    print(f"{len(wanted) - count} IDs from {args.filter} not found in {args.input}")
print(f"Saved {count} records from {args.input} to {args.output}")



