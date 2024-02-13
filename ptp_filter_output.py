import csv
import argparse

# Argument parser
parser = argparse.ArgumentParser(description="Search GenBank file, retrieve gene sequences and save as fasta.")
parser.add_argument("-i", "--input", type=str, help="mPTP output txt file")
parser.add_argument("-m", "--meta", type=str, help="metadata file (supermatrix_count.py output)")
parser.add_argument("-k", "--keep", type=str, help="file with list of taxa to keep (eg constraint, outgroup)")
parser.add_argument("-o", "--output", type=str, help="output file with filtered fasta IDs")

args = parser.parse_args()

# Save gene count for each taxon in dict
count = {}
with open(args.meta) as file:
    metadata = csv.reader(file)
    for row in metadata:
        try:
            count[row[0]] = int(row[2])
        except ValueError:
            continue

keep = []
if args.keep:
    file = open(args.keep)
    lines = file.readlines()
    for line in lines:
        keep.append(line.strip)

file = open(args.input)

# Save mPTP delimited species lists
species_lists = []
temp = []
x = 0
lines = file.readlines()
for line in lines:
    if 'Number of delimited species' in line:
        print(line)
    x += 1
    if x < 10:
        continue
    if 'Species ' in line:
        continue
    line = line.strip()
    if line != '':
        temp.append(line)
    else:
        species_lists.append(temp)
        temp = []
species_lists.append(temp)

# Get taxon with most genes in each species list
chosen = []
for species in species_lists:
    nt = 0   # nucleotide count
    ch = ''  # chosen taxon
    lab = 0
    for s in species:
        if s in keep:
            chosen.append(s)
            lab = 1
        else:
            try:
                if count[s] > nt:
                    nt = count[s]
                    ch = s
            except KeyError:
                print(f'{s} is not in metadata')

    if lab == 0:
        chosen.append(ch)

output = open(args.output, 'w')
for m in chosen:
    if m != '':
        output.write(f'{m}\n')

print(f'{len(chosen)} taxa saved to {args.output}')
