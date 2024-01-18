
#python3 get_concat_recs.py -e mixedupvoyage@gmail.com -t Megadytes

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Search GenBank file, retrieve gene sequences and save as fasta.")
parser.add_argument("-i", "--input", type=str, help="input acc list")
parser.add_argument('-g', '--gb_file', type=str, help="Input genbank format file")
parser.add_argument("-o", "--output", type=str, help="output txid list file")


args = parser.parse_args()

accs = []
file = open('test.txt')
lines = file.readlines()
for line in lines:
    acc = line.strip()
    accs.append(acc)
print(f'{len(accs)} IDs found')


# Search through GBIDs
x = 0  # Count taxids
txids = []
with open('test.gb') as file:
    record = SeqIO.parse(file, "gb")
    for rec in record:
        if rec.name in accs:
            try:
                db_xref = rec.features[0].qualifiers["db_xref"]
            except KeyError:
                continue
            for ref in db_xref:
                if "taxon" in ref:  # Get NCBI taxon, rather than BOLD cross ref
                    txid = "".join(filter(str.isdigit, ref))  # Extract numbers from NCBI taxon value
            txids.append(txid)

with open(args.output, 'w') as file:
    for txid in txids:
        file.write(f'{txid}\n')