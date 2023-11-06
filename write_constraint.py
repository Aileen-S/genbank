from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Filter fasta file from list of IDs. Separate multiple files with commas.")
parser.add_argument("-i", "--ingroup", type=str, help="Ingroup fasta")
parser.add_argument("-o", "--outgroup", type=str, help="Outgroup fasta")
args = parser.parse_args()

ins = []
ingroup = args.ingroup.split(',')
for file in ingroup:
    for taxon in SeqIO.parse(file, "fasta"):
        ins.append(taxon.id)
        #fast.write(f'>{taxon.id}\n1\n')

outs = []
outgroup = args.outgroup.split(',')
for file in outgroup:
    for taxon in SeqIO.parse(file, "fasta"):
        outs.append(taxon.id)
        #fast.write(f'>{taxon.id}\n0\n')

with open('fasttree_constraint.txt', 'w') as fast:
    for i in ins:
        fast.write(f'{i}\n1\n')
    for o in outs:
        fast.write(f'{o}\n0\n')

with open('raxml_constraint.txt', 'w') as rax:
    rax.write(f'(({",".join(outs)}),({",".join(ins)}));')

with open('outgroup.txt', 'w') as file:
    for o in outs:
        file.write(f'{o},')
