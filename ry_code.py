from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="RY coding")
parser.add_argument("-i", "--input", type=str, help="Input fasta")
parser.add_argument("-o", "--output", type=str, help="Output fasta")
args = parser.parse_args()

with open(args.output, "w") as output:
    records = SeqIO.parse(args.input, "fasta")
    for rec in records:
        seq = rec.seq.upper()
        seq = seq.replace('A', 'R').replace('G', 'R').replace('C', 'Y').replace('T', 'Y')
        output.write(f'>{rec.id}\n{seq}\n')
