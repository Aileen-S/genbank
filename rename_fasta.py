from Bio import SeqIO
import argparse
import csv

parser = argparse.ArgumentParser(description="Rename sequences in fasta file from CSV.")
parser.add_argument("-i", "--input", type=str, help="Input fasta")
parser.add_argument("-n", "--names", type=str, help="List of new names (with old name as part of name")
parser.add_argument("-o", "--output", type=str, help="Output file")
args = parser.parse_args()

new = []
with open(args.names) as file:
    lines = file.readlines()
    for line in lines:
        line = line.strip()
        new.append(line)

records = SeqIO.parse(args.input, "fasta")
output = open(args.output, 'w')
for rec in records:
    if ';frame=' in rec.id:
        recid, frame = rec.id.split(';')
        for n in new:
            if recid in n:
                rec.id = f'{n};{frame}'
        output.write(f'>{rec.id}\n{rec.seq}\n')
    else:
        for n in new:
            if rec.id in n:
                rec.id = n
        output.write(f'>{rec.id}\n{rec.seq}\n')

print(f'Saved renamed fasta to {args.output}')