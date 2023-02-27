
import argparse

parser = argparse.ArgumentParser(description="Get lab mitogenome genbanks from ID list")
parser.add_argument("-i", "--input", type=str, help="Input ID list file")
parser.add_argument("-o", "--output", type=str, help="Name output file: default is <inputfile>.out")
parser.add_argument('-v', '--version', type=str, help='Database version, eg: gbmaster_2022-06-27')

args = parser.parse_args()

# Use file with list of dbids to combine individual genbank records into one file.
file = open(args.input)
if args.output:
    output = open(args.output, "w")
else:
    output = open(f'{args.input}.out', 'w')

x = 0
lines = file.readlines()
print(f'Searching /mbl/share/workspaces/groups/voglerlab/MMGdatabase/{args.version}/ for IDs in {args.input}')
for line in lines:
    line = line.strip()
    try:
        name = f'/mbl/share/workspaces/groups/voglerlab/MMGdatabase/{args.version}/{line}.gb'
        x += 1
    except FileNotFoundError:
        print(f'No record found for {line}')
        continue
    record = open(name)
    record = record.read()
    output.write(record)

print(f'{x} records saved to {output}')
