import argparse

parser = argparse.ArgumentParser(description="Create RAxML partition files from IQTREE best_scheme.nex")
parser.add_argument("-i", "--input", type=str, help="in file")
#parser.add_argument('-t', '--sequence_type', choices=['aa', 'nt'], help="Is input file from amino acid or nucleotide supermatrix.")

args = parser.parse_args()

mods = {'JC': ['JC69'],
          'K80': ['K2P'],
          'HKY': ['HKY85'],
          'TrN': ['TN', 'TN93'],
          'TrNef': ['TNe'],
          'TPM1': ['K81', 'K3P'],
          'TPM1uf': ['K81u'],
          'TPM2uf': ['TPM2u'],
          'TPM3uf': ['TPM3u'],
          'TIM1uf': ['TIM'],
          'TIM1': ['TIMe'],
          'TIM2uf': ['TIM2'],
          'TIM2': ['TIMe'],
          'TIM3uf': ['TIM3'],
          'TIM3': ['TIMe'],
          'TVMuf': ['TVM'],
          'TVM': ['TVMe'],}

file = open(args.input)
lines = file.readlines()

models = {}

x = 0
y = 0
output = open(f'{args.input}_RAxML.txt', 'w')
for line in lines:
    line = line.strip(';\n')
    if 'charset' in line:
        y += 1
        split = line.split(' = ')
        name = split[0].replace('  charset ', '')
        location = split[1].replace('  ', ',')
        models[y] = [name, location]
    if ': ' in line:
        x += 1
        split = line.split(': ')
        model = split[0].strip('    ').split('+', 1)
        for key, value in mods.items():
            for v in value:
                if v == model[0]:
                    model[0] = model[0].replace(v, key)
        model = '+'.join(model)
        models[x].append(model)

for v in models.values():
    output.write(f'{v[2]}, {v[0]}={v[1]}\n')

