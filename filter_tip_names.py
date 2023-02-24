# Get ID list to use as filter
#ids = []
#file = open('test.txt')
#lines = file.readlines()
#for line in lines:
#    line = line.strip()
#    ids.append(line)
#print(ids)

import argparse

parser = argparse.ArgumentParser(description="Get list of accessions or taxon IDs from fasta.")
parser.add_argument("-i", "--input", type=str, help="Input file to be filtered")
parser.add_argument("-f", "--filter", type=str, help="File with list of IDs to search for in input file")
args = parser.parse_args()

infile = open(args.input)
filter = open(args.filter)
output = open(f'{args.input}.out', 'w')

ids = []
lines = filter.readlines()
for line in lines:
    line = line.strip()
    ids.append(line)

lines = infile.readlines()
for line in lines:
    line = line.strip()
    for i in ids:
        if i in line:
            output.write(f'{line}\n')