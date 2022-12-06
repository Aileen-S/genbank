import argparse

parser = argparse.ArgumentParser(description="Create RAxML partition files from IQTREE best_scheme.nex")
parser.add_argument("-i", "--input", type=str, help="in file")
#parser.add_argument('-t', '--sequence_type', choices=['aa', 'nt'], help="Is input file from amino acid or nucleotide supermatrix.")

args = parser.parse_args()

file = open(args.input)
lines = file.readlines()

models = {}

x = 0
y = 0
output = open(f'{args.input}_RAxML.txt', 'w')
for line in lines:
    line = line.strip()
    if 'charset' in line:
        y += 1
        split = line.split(' = ')
        name = split[0].replace('charset ', '')
        location = split[1].replace('  ', ', ')
        models[y] = [name, location]
    if ': ' in line:
        x += 1
        split = line.split(': ')
        models[x].append(split[0])

for v in models.values():
    output.write(f'{v[2]}, {v[0]}={v[1]}\n')
    #start = split[0].split('/')
    #start[0] = 'DNA'
    #output.write(f'{start[0]}, {start[1]} = {split[1]}\n')
