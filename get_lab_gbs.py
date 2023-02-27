
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

lines = file.readlines()
for line in lines:
    line = line.strip()
    name = f'/mbl/share/workspaces/groups/voglerlab/MMGdatabase/{args.version}/{line}.gb'
    record = open(name)
    record = record.read()
    output.write(record)
