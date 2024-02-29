from Bio import Entrez
import argparse

# Argument parser
parser = argparse.ArgumentParser(description="Download GenBank files for requested taxon.")
parser.add_argument("-t", "--taxon", type=str, help="Taxon of interest")
parser.add_argument("-e", "--email", type=str, help="Your email registered with NCBI")
parser.add_argument("-o", "--output", type=str, help="Output file name")

args = parser.parse_args()         # Process input args from command line

Entrez.email = args.email

output = open(args.output, 'w')

# Get initial count of records
handle = Entrez.esearch(db="nucleotide", term=args.taxon, retmax=0)
record = Entrez.read(handle)
count = int(record["Count"])
print(str(count) + " records found")

# Get list of GenBank ID numbers for records
gbids = []
x = 0
y = 0
for start in range(0, count, 1000):
    print(f"Getting GenBank IDs for taxon IDs {x + 1} to {x + 1000}" if (x + 1000) < count else
          f"Getting GenBank IDs for taxon IDs {x + 1} to {count}")
    x += 1000
    y += 1
    # Search and get GB IDs
    handle = Entrez.esearch(db="nucleotide", term=args.taxon, retstart=start, retmax=1000)
    record = Entrez.read(handle)
    gbids.append(','.join(record['IdList']))

# Use GBIDs to download GenBank records
z = 0
with open(args.output, 'a') as file:
    for gbid_str in gbids:
        print(f"Downloading GenBank records for taxon IDs {z + 1} to {z + 1000}" if (z + 1000) < count else
              f"Downloading GenBank records for taxon IDs {z + 1} to {count}")
        z += 1000
        #if z <= 24000:
        #    continue
        handle = Entrez.efetch(db="nucleotide", id=gbid_str, rettype="gb", retmode="text")
        output.write(handle.read())
print(f'Output saved to {args.output}')

