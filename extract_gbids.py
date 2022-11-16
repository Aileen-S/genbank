import argparse

parser = argparse.ArgumentParser(description="Get list of accessions from fasta.")
parser.add_argument("-i", "--input", type=str, help="in file")
#parser.add_argument("-o", "--output", type=str, help="out file")

args = parser.parse_args()

file = open(args.input)
output = open(f'{args.input}.out', 'w')
lines = file.readlines()
for line in lines:
    line.strip()
    #if ">" not in line:
        #continue
    name = line.split("_")
    if "NC_" in line:
        acc = f"{name[1]}_{name[2]}"
    else:
        acc = name[1]
    output.write(f'{acc}\n')


