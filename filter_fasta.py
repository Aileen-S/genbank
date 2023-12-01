from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Filter fasta file from list of IDs")
parser.add_argument("-i", "--input", type=str, help="Input fasta to be filtered")
parser.add_argument("-f", "--found", type=str, help="Output file for sequences in search filter")
parser.add_argument("-n", "--notfound", type=str, help="Output file for input sequences not present in search filter")
parser.add_argument("-s", "--search", type=str, help="File with IDs to search for. Specify file type with -t.")
parser.add_argument('-t', '--type', type=str, choices=['list', 'fasta', 'hmmer'], help="Search file type:\n"
                                                                                         "list: file with list of complete or partial fasta IDs\n"
                                                                                         "fasta: fasta file\n"
                                                                                          "hmmer: file with HMMER output from --tblout option")
args = parser.parse_args()         # Process input args from command line
#args = argparse.Namespace(input='ATP6.fasta', filter='test.txt')

wanted = []
with open(args.search) as infile:
    if args.type == 'list':
        lines = infile.readlines()
        for line in lines:
            wanted.append(line.strip())

    elif args.type == 'fasta':
        recs = SeqIO.parse(args.search, "fasta")
        for rec in recs:
            wanted.append(rec.id)

    elif args.type == 'hmmer':
        lines = infile.readlines()
        for line in lines:
            if '#' in line:
                continue
            line = line.split(' ')
            wanted.append(line[0])
print(f"Found {len(wanted)} unique identifiers in {args.search}")

print(f"Searching {args.input}")

if args.found:
    records = (r for r in SeqIO.parse(args.input, "fasta") if r.id in wanted)
    count = SeqIO.write(records, args.found, "fasta")
    print(f"{count} of {len(wanted)} records from {args.search} found in {args.input}")
    print(f"Saved {count} records from {args.input} to {args.found}")

if args.notfound:
    records = (r for r in SeqIO.parse(args.input, "fasta") if r.id not in wanted)
    count = SeqIO.write(records, args.notfound, "fasta")
    print(f"{count} records from {args.input} not present in {args.search}")
    print(f"Saved {count} records from {args.input} to {args.notfound}")



