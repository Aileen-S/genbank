from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Filter fasta file from list of IDs. Choose either -p, -s, -h or -f flag.")
parser.add_argument("-i", "--ingroup", type=str, help="Ingroup fasta")
parser.add_argument("-o", "--outgroup", type=str, help="Outgroup fasta")
parser.add_argument("-f", "--outfile", type=str, help="Output file")

ins = []
args = parser.parse_args()
for taxon in SeqIO.parse(args.ingroup, "fasta"):
    ins.append(taxon.id)

outs = []
for taxon in SeqIO.parse(args.outgroup, "fasta"):
    outs.append(taxon.id)

with open(args.outfile, 'w') as file:
    file.write(f'(({",".join(outs)}),({",".join(ins)}));')