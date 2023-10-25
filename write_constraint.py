from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Filter fasta file from list of IDs. Choose either -p, -s, -h or -f flag.")
parser.add_argument("-i", "--ingroup", type=str, help="Ingroup fasta")
parser.add_argument("-o", "--outgroup", type=str, help="Outgroup fasta")
args = parser.parse_args()

fast = open("fasttree_constraint.txt", 'w')
rax = open('raxml_constraint.txt', 'w')
ins = []

for taxon in SeqIO.parse(args.ingroup, "fasta"):
    ins.append(taxon.id)
    fast.write(f'>{taxon.id}\n1\n')

outs = []
for taxon in SeqIO.parse(args.outgroup, "fasta"):
    outs.append(taxon.id)
    fast.write(f'>{taxon.id}\n0\n')

with open('raxml_constraint.txt', 'w') as rax:
    rax.write(f'(({",".join(outs)}),({",".join(ins)}));')

