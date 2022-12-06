import argparse

parser = argparse.ArgumentParser(description="Get list of accessions or taxon IDs from fasta.")
parser.add_argument("-i", "--input", type=str, help="in file")
parser.add_argument("-r", "--ref", choices=['txid', 'gbid'], help="Choose GBID or TXID")

#parser.add_argument("-o", "--output", type=str, help="out file")

args = parser.parse_args()

file = open(args.input)
output = open(f'{args.input}.out', 'w')
lines = file.readlines()
for line in lines:
    line.strip()
    if line.startswith('>'):
        name = line.split("_")
        if args.ref == 'gbid':
            if "NC_" in line:
                gbid = f"{name[1]}_{name[2]}"
            else:
                gbid = name[1]
            output.write(f'{gbid}\n')
        if args.ref == 'txid':
            txid = name[0].replace('>', '')
            output.write(f'{txid}\n')
    else:
        continue


