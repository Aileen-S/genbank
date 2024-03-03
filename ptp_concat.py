
import csv
import argparse

# Argument parser
parser = argparse.ArgumentParser(description="Process PTP output file.")
parser.add_argument("-i", "--input", type=str, help="mPTP output txt file")
parser.add_argument("-m", "--meta", type=str, help="metadata file (supermatrix_count.py output)")
parser.add_argument("-k", "--keep", type=str, help="file with list of taxa to keep (eg constraint, outgroup)")
parser.add_argument("-o", "--output", type=str, help="output file with filtered fasta IDs")

args = parser.parse_args()

# Save gene count for each taxon in dict
count = {}
meta = {}
with open(args.meta) as file:
    metadata = csv.reader(file)
    for row in metadata:
        try:
            count[row[0]] = int(row[2])
            meta[row[0]] = row
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
ptp_meta = {}
species_lists = []
temp = []
x = 0
y = 0
lines = file.readlines()
for line in lines:
    y += 1
    if 'Number of delimited species' in line:
        print(line)
    x += 1
    if y < 8:
        continue
    if 'Species ' in line:
        spec_no = line.strip().replace(':', '')
        spec_no = spec_no.split(' ')[1]
        spec_no = str(spec_no).zfill(4)
        spec_no = f'Taxon_{spec_no}'
        continue
    line = line.strip()
    if line != '':
        if spec_no in ptp_meta:
            ptp_meta[spec_no].append(meta[line])
        else:
            ptp_meta[spec_no] = [meta[line]]
        temp.append(line)
    else:
        species_lists.append(temp)
        temp = []
species_lists.append(temp)


# Write PTP metadata
with open('ptp_metadata.csv', 'w') as file:
    file.write('PTP Species,Taxon,Total Genes,Total Nucleotides,12S,16S,18S,28S,AK,ATP6,ATP8,CAD,COX1a,COX1b,COX2,COX3,CYTB,EF1A,H3,ND1,ND2,ND3,ND4,ND4L,ND5,ND6,RNApol,Wg\n')
    writer = csv.writer(file)
    for k, values in ptp_meta.items():
        x += 1
        for v in values:
            row = [k] + v
            writer.writerow(row)

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
