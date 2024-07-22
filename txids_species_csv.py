from Bio import Entrez
import csv
import argparse
import time

# Get species names from TXID list
def fetch_species_names(taxon_ids):
    taxon_data = []
    
    for taxon_id in taxon_ids:
        for attempt in range(1, 10):
            try:
                handle = Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml")
                records = Entrez.read(handle)
                handle.close()
            except Entrez.HTTPError:
                print("HTTP error fetching records; retry in 20 seconds")
                time.sleep(20)
            else:
                break
        else: 
            'Failed to retrieve records. Try again later.'
        if records:
            species_binomial = records[0]["ScientificName"]
            rank = records[0]["Rank"]
            taxon_data.append((taxon_id, species_binomial, rank))
        else:
            taxon_data.append((taxon_id, "Not Found"))
    
    return taxon_data


parser = argparse.ArgumentParser(description="Get csv with txids and binomials")
parser.add_argument("-e", "--email", type=str, help="Your email registered with NCBI")
parser.add_argument('-i', '--input', type=str, help='file with NCBI taxon ID list')
parser.add_argument("-o", "--output", type=str, help="Output csv")

args = parser.parse_args()

# Set your email address
Entrez.email = args.input


taxids = []
file = open(args.input)
lines = file.readlines()
for line in lines:
    taxid = line.strip()
    taxids.append(taxid)
print(f'{len(taxids)} IDs found in {args.input}')

taxids = ['3029442', '3031963', '3039960', '3046642', '3048115', '3053695', '3053696', '3078427']

# Fetch species names
taxon_data = fetch_species_names(taxids)
print(taxon_data)

# Write to CSV
with open(args.output, "w") as output:
    csvwriter = csv.writer(output)
    csvwriter.writerow(["ncbi_taxid", "species", "rank"])
    csvwriter.writerows(taxon_data)

print(f"Printed output to {args.output}")
