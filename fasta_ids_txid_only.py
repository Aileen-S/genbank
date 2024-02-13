import argparse

parser = argparse.ArgumentParser(description="Stip fasta to TXID if available, other ID if not.")
parser.add_argument("-i", "--input", type=str, help="in file")
parser.add_argument("-o", "--output", type=str, help="out file")
args = parser.parse_args()

exceptions = ['con_', 'lab_', 'outgroup_']

file = open(args.input)
output = open(args.output, 'w')
lines = file.readlines()
for line in lines:
    try:
        if line.startswith('>'):
            line = line.replace('>', '')
            name = line.split("_")
            x = 0
            for n in name:
                if not any(char.isalpha() for char in n):
                    if len(n) > 4:
                        line = f'{n}\n'
                        x =+ 1
                        continue
            if x == 0:
                for x in exceptions:
                    if x in line:
                        print(line)
                        line = line.replace(x, '')
                line = line.split('_')
                print(line)
                line = f'{line[0]}\n'
            line = f'>{line}'
    except IndexError:
        pass
    output.write(line)


