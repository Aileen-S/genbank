import csv
import argparse

# Argument parser
parser = argparse.ArgumentParser(description="Search GenBank file, retrieve gene sequences and save as fasta.")
parser.add_argument("-i", "--input", type=str, help="mPTP output txt file")
parser.add_argument("-m", "--meta", type=str, help="metadata file")
parser.add_argument("-o", "--output", type=str, help="output file with filtered fasta IDs")

args = parser.parse_args()

# Save gene count for each taxon in dict
count = {}
with open(args.meta) as file:
    metadata = csv.reader(file)
    for row in metadata:
        count[row[0]]: row[5]
print(count)

file = open(args.input)

# Save mPTP delimited species lists
species_lists = []
temp = []
lines = file.readlines()
for line in lines:
    line = line.strip()
    if line != '':
        temp.append(line)
    else:
        species_lists.append(temp)
        temp = []
species_lists.append(temp)

# Get taxon with most genes in each species list
most = []
for species in species_lists:
    c = 0
    spec = ''
    for s in species:
        if 'Species' in s:
            continue
        if len(s) <= 15:
            most.append(s)
        else:
            if s in count:
                if count[s] > 0:
                    spec = s
        most.append(spec)

output = open(args.output, 'w')
for m in most:
    if m != '':
        output.write(f'{m}\n')