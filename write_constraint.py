from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Filter fasta file from list of IDs. Choose either -p, -s, -h or -f flag.")
parser.add_argument("-i", "--ingroup", type=str, help="Ingroup fasta")
parser.add_argument("-o", "--outgroup", type=str, help="Outgroup fasta")
parser.add_argument("-f", "--outfile", type=str, help="Output file")

ins = []
args = parser.parse_args()
ingroup = SeqIO.parse(args.ingroup, "fasta")
for i in ingroup:
    ins.append(i.id)

outs = []
args = parser.parse_args()
outgroup = SeqIO.parse(args.outgroup, "fasta")
for i in outgroup:
    outs.append(i.id)

with open(args.outfile, 'w') as file:
    file.write(f'(({",".join(outgroup)}),({",".join(ingroup)}));')