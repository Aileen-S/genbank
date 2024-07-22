from Bio import Entrez
import csv
import argparse


def fetch_species_names(taxon_ids):
    taxon_data = []
    
    for taxon_id in taxon_ids:
        handle = Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        if records:
            species_binomial = records[0]["ScientificName"]
            taxon_data.append((taxon_id, species_binomial))
        else:
            taxon_data.append((taxon_id, "Not Found"))
    
    return taxon_data


parser = argparse.ArgumentParser(description="Get csv with txids and binomials")
parser.add_argument("-e", "--email", type=str, help="Your email registered with NCBI")
parser.add_argument('-i', '--input', action=str, help='file with NCBI taxon ID list')
parser.add_argument("-o", "--output", type=str, help="Output csv")

args = parser.parse_args()

# Set your email address
Entrez.email = args.input


taxids = []
file = open(args.file)
lines = file.readlines()
for line in lines:
    taxid = line.strip()
    taxids.append(taxid)
print(f'{len(taxids)} IDs found in {args.inout}')


# Fetch species names
taxon_data = fetch_species_names(taxids)

# Write to CSV
with open(args.output, "w") as output:
    csvwriter = csv.writer(output)
    csvwriter.writerow(["ncbi_taxid", "species"])
    csvwriter.writerows(taxon_data)

print(f"Printed output to {args.output}")
