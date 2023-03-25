from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Filter fasta file from list of IDs. Choose either -p, -s, -h or -f flag.")
parser.add_argument("-i", "--input", type=str, help="Input fasta to be filtered")
parser.add_argument("-o", "--output", type=str, help="Output file")
parser.add_argument("-d", "--duds", action='store_true', help="Write separate fasta with records not found with filter")
parser.add_argument("-do", "--dudsout", action='store_true', help="Output file for rejected sequences")

# Filter options
parser.add_argument('-s', '--string', type=str, help='Single string to search for in sequence identifiers')
parser.add_argument('-p', '--partial_file', type=str, help='File containing TXIDs to search for in sequence identifiers')
parser.add_argument("-f", "--filter", type=str, help="File containing list of sequence identifiers to extract from input file")
parser.add_argument("-hm", "--hmmer", type=str, help="File containing HMM output from --tblout flag")

args = parser.parse_args()

#args = argparse.Namespace(input='ATP6.fasta', filter='test.txt')

if args.partial_file:
    with open(args.partial_file) as infile:
        ids = set(line.rstrip("\n").split(None, 1)[0] for line in infile)
        print(f"Found {len(ids)} unique identifiers in {args.partial_file}")
        wanted = []
        for r in SeqIO.parse(args.input, "fasta"):
            for i in ids:
                if i in r.id:
                    wanted.append(r.id)

elif args.string:
    wanted = []
    ids = [args.string]
    for r in SeqIO.parse(args.input, "fasta"):
        for i in ids:
            if i in r.id:
                wanted.append(r.id)

elif args.hmmer:
    wanted = []
    with open(args.hmmer) as infile:
        lines = infile.readlines()
        for line in lines:
            if '#' in line:
                continue
            line = line.split(' ')
            wanted.append(line[0])
            print(f"Found {len(wanted)} unique identifiers in {args.hmmer}")



else:
    with open(args.filter) as infile:
        wanted = set(line.rstrip("\n") for line in infile)
        print(f"Found {len(wanted)} unique identifiers in {args.filter}")
print(f"Searching {args.input}")

records = (r for r in SeqIO.parse(args.input, "fasta") if r.id in wanted)
count = SeqIO.write(records, args.output, "fasta")

if count < len(wanted):
    print(f"{len(wanted) - count} IDs from {args.filter} not found in {args.input}")
print(f"Saved {count} records from {args.input} to {args.output}")

if args.duds:
    records = (r for r in SeqIO.parse(args.input, "fasta") if r.id not in wanted)
    count = SeqIO.write(records, args.dudsout, "fasta")
    print(f"Saved {count} records from {args.input} to {args.duds}")



