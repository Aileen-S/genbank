import argparse

parser = argparse.ArgumentParser(description="Filter fasta file from list of IDs. Choose either -p, -s, -h or -f flag.")
parser.add_argument("-i", "--input", type=str, help="HMM output file from --tblout flag")
parser.add_argument("-o", "--output", type=str, help="Output file with list of IDs")
args = parser.parse_args()

#args = argparse.Namespace(input='ATP6.fasta', filter='test.txt')

output = open(args.output, 'w')
x = 0
with open(args.input) as infile:
    lines = infile.readlines()
    for line in lines:
        if '#' in line:
            continue
        line = line.split('.', 1)
        output.write(f'{line[0]}\n')
        x += 1
print(f'{x} accessions saved to {args.output}')