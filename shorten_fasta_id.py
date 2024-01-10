from Bio import SeqIO
import argparse
import csv

parser = argparse.ArgumentParser(description="Rename sequences in fasta file from CSV.")
parser.add_argument("-i", "--input", type=str, help="Input fasta")
parser.add_argument("-o", "--output", type=str, help="Output file")
parser.add_argument("-n", "--number", type=int, help="Strip fasta ID to this many components")

args = parser.parse_args()


records = SeqIO.parse(args.input, "fasta")
output = open(args.output, 'w')
for rec in records:
    if ';frame=' in rec.id:
        recid, frame = rec.id.split(';')
        name = recid.split('_')
        rec.id = f'{"_".join(name[0:args.number])};{frame}'
        output.write(f'>{rec.id}\n{rec.seq}\n')

    else:
        name = rec.id.split('_')
        rec.id = f'{"_".join(name[0:3])}'
        output.write(f'>{rec.id}\n{rec.seq}\n')

print(f'Saved renamed fasta to {args.output}')