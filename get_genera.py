

#python3 get_genera.py -e mixedupvoyage@gmail.com -t Eretes -o test.out

import argparse
from Bio import Entrez
from Bio import SeqIO

# Argument parser
parser = argparse.ArgumentParser(description="Get list of genera represented on GenBank. Input either higher taxonomy Latin name or taxon ID.")
parser.add_argument("-t", "--taxon", type=str, help="Taxon of interest - Linean name or NCBI higher taxon ID")
parser.add_argument("-o", "--output", type=str, help="Output file")
parser.add_argument("-e", "--email", type=str, help="Your email registered with NCBI")
args = parser.parse_args()         # Process input args from command line
Entrez.email = args.email

if args.taxon.isdigit():
    print('yes')
    taxon = f'txid{args.taxon}[Organism:exp]'
else:
    taxon = args.taxon
print(taxon)

# Get initial count of records
handle = Entrez.esearch(db="nucleotide", term=taxon, retmax=0)
record = Entrez.read(handle)
count = int(record["Count"])
print(str(count) + " records found")

genera = set()
#gbids = []
x = 0
y = 0
for start in range(0, count, 10000):
    print(f"Searching GenBank records {x + 1} to {x + 10000}" if (x + 10000) < count else
          f"Searching GenBank records {x + 1} to {count}")
    x += 10000
    y += 1
    # Search and get GB IDs
    handle = Entrez.esearch(db="nucleotide", term=taxon, retstart=start, retmax=10000)
    record = Entrez.read(handle)
    gbids = ','.join(record['IdList'])
    handle = Entrez.efetch(db="nucleotide", id=gbids, rettype="gb", retmode="text")  # Get GenBanks
    record = SeqIO.parse(handle, "gb")
    for rec in record:
        gen_spec = rec.annotations["organism"].split(' ')
        genus = gen_spec[0]
        genera.add(genus)

print(f'{len(genera)} genus found' if len(genera) == 1 else
      f'{len(genera)} genera found')


output = open(args.output, 'w')
for genus in genera:
    output.write(f'{genus}\n')
print(f'Genus save to {args.output}' if len(genera) == 1 else
      f'Genera save to {args.output}')