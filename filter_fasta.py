from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Filter fasta file from list of IDs")
parser.add_argument("-i", "--input", type=str, help="Input fasta to be filtered")
parser.add_argument("-f", "--found", type=str, help="Output file for sequences in search filter")
parser.add_argument("-n", "--notfound", type=str, help="Output file for input sequences not present in search filter")
parser.add_argument("-s", "--search", type=str, help="File with IDs to search for. Specify file type with -t.")
parser.add_argument("-p", "--partial", action='store_true', help="Search file contains only part of input fasta ids (eg without frame tags or taxonomy)")

parser.add_argument('-t', '--type', type=str, choices=['list', 'fasta', 'hmmer'],
                    help="Search file type:\n"
                         "list: file with list of complete or partial fasta IDs\n"
                         "fasta: fasta file\n"
                         "hmmer: file with HMMER output from --tblout option\n")

args = parser.parse_args()         # Process input args from command line

filter = []
with open(args.search) as infile:
    if args.type == 'list':
        lines = infile.readlines()
        for line in lines:
            filter.append(line.strip())

    elif args.type == 'fasta':
        recs = SeqIO.parse(args.search, "fasta")
        for rec in recs:
            filter.append(rec.id)

    elif args.type == 'hmmer':
        lines = infile.readlines()
        for line in lines:
            if '#' in line:
                continue
            line = line.split(' ')
            filter.append(line[0])

len_filter = len(filter)
print(f'{len_filter} IDs in {args.search}')

records = SeqIO.parse(args.input, "fasta")
print(f'{len(list(records))} records in {args.input}')

if args.partial:
    filter_ids = []
    records = SeqIO.parse(args.input, "fasta")
    x = 0
    for rec in records:
        for w in filter:
            if w in rec.id:
                filter_ids.append(rec.id)
    filter = filter_ids

if args.found:
    records = (r for r in SeqIO.parse(args.input, "fasta") if r.id in filter)
    count = SeqIO.write(records, args.found, "fasta")
    print(f"{count} of {len_filter} records from {args.search} found in {args.input}")
    print(f"{count} records from {args.input} saved to {args.found}")

if args.notfound:
    records = (r for r in SeqIO.parse(args.input, "fasta") if r.id not in filter)
    count = SeqIO.write(records, args.notfound, "fasta")
    print(f"{count} records from {args.input} not present in {args.search}")
    print(f"{count} records from {args.input} saved to {args.notfound}")

print('\n')
