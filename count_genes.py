
#python3 count_genes.py -e mixedupvoyage@gmail.com -t Eretes

import argparse
import csv
from Bio import Entrez
from Bio import SeqIO


# Function definitions

def get_feat_name(feat):
    featname = "unknown"
    nametags = ['gene', 'product', 'label', 'standard_name']  # Search these four keys for gene name
    if any(t in feat.qualifiers.keys() for t in nametags):
        for t in nametags:
            if t in feat.qualifiers.keys():
                featname = feat.qualifiers[t][0].upper()
                break
    return featname


def set_feat_name(feat, name):
    nametags = ['gene', 'product', 'label', 'standard_name']
    if any(t in feat.qualifiers.keys() for t in nametags):
        for t in nametags:
            if t in feat.qualifiers.keys() :
                feat.qualifiers[t][0] = name
    return feat


def search_nuc(term, summaries=False, chunk=10000):
    # Get initial count of responses
    searchhand = Entrez.esearch(db="nucleotide", term=term, retmax=0)
    searchrec = Entrez.read(searchhand)
    count = int(searchrec["Count"])
    print(str(count) + " records found")
    # Yield
    for start in range(0, count, chunk):
        # Search and get GB IDs
        searchhand = Entrez.esearch(db="nucleotide", term=term, retstart=start, retmax=chunk)
        searchrec = Entrez.read(searchhand)
        gbids = searchrec['IdList']
        # Yield only GBIDs if no summaries desired
        if not summaries:
            yield gbids
        else:
            # Retrieve summaries and yield both otherwise
            sumhand = Entrez.esummary(db="nucleotide", id=','.join(gbids))
            sumrec = Entrez.read(sumhand)
            yield gbids, sumrec


# Argument parser
parser = argparse.ArgumentParser(description="Search GenBank, count records for each gene name")
parser.add_argument("-t", "--taxon", type=str, help="Taxon of interest")
parser.add_argument("-e", "--email", type=str, help="Your email registered with NCBI")


# Start the actual script
args = parser.parse_args()         # Process input args from command line
#args = argparse.Namespace(taxon='Amphizoidae', mpc=True, email='aileen.scott@nhm.ac.uk', nuclear=False) # This is how I step through the script interactively
#Namespace(taxon='Eretes', mpc=True, email='aileen.scott@nhm.ac.uk', nuclear=False)


Entrez.email = args.email

# Generate search term to get all sequences in the search taxonomy
# - if -n option not used, then include "mitochondrial" in search term.
basesearch = f"(\"{args.taxon}\"[Organism] OR \"{args.taxon}\"[All Fields])"

# Retrieve all taxids that represent tips of the NCBI Taxonomy tree
# Make the search generator
searchgen = search_nuc(term=basesearch, summaries=True, chunk=5000)

taxids = set()
i = 0
for gbids, summaries in searchgen:
    # gbids, summaries = next(searchgen)
    i += 1
    taxa = set(int(s['TaxId']) for s in summaries)
    taxids.update(taxa)
    print(f"iteration={i}, returns={len(gbids)}, first gbid={gbids[0]}, first summary accession={summaries[0]['Caption']}, taxids in this iteration={len(taxa)}, total taxids={len(taxids)}")

# Some of these will be subspecies.
# You need to search them in NCBI Taxonomy to weed out the subspecies and generate a list of latin biomials.
# Then iterate through each of these binomials (not taxids as initially thought) to download the sequences etc


print(f"{len(taxids)} unique taxon IDs saved")
print("Searching GenBank")
print("Downloading GenBank records for taxon IDs 0 to 100" if len(taxids) > 100 else
      f"Downloading GenBank records for taxon IDs 0 to {len(taxids)}")

# Set accepted genes and minimum sequence lengths
#min = {"ATP6": 500, "ATP8": 100, "COX1": 500, "COX2": 500, "COX3": 500, "CYTB": 500, "ND1": 500, "ND2": 500, "ND3": 300, "ND4": 500, "ND4L": 200, "ND5": 500, "ND6": 400}

x = 0  # Count taxids
y = 0  # Count records saved
genes = {}
for tax in taxids:
    y += 1
    if y % 100 == 0:
        print(f"Downloading GenBank records for taxon IDs {y+1} to {y+100}" if (y+100) < len(taxids) else
              f"Downloading GenBank records for taxon IDs {y+1} to {len(taxids)}")
    handle = Entrez.esearch(db="nucleotide", term=f"txid{tax}")       # Search for all records for each taxon id
    record = Entrez.read(handle)
    accs   = record["IdList"]                                         # Get accessions
    accstr = ",".join(accs)                                           # Join into string for efetch
    handle = Entrez.efetch(db="nucleotide", id=accstr, rettype="gb", retmode="text")  # Get GenBanks
    record = SeqIO.parse(handle, "gb")
    for rec in record:
        if args.taxon not in rec.annotations["taxonomy"]:
            continue
        if len(rec.annotations["taxonomy"]) >= 14:
            subfamily = rec.annotations["taxonomy"][13]
        else:
            subfamily = undefined
        for feature in rec.features:
            type = feature.type
            if type not in ('CDS', 'rRNA', 'mRNA'):
                continue  # skip the rest of the current iteration of this loop
            name = get_feat_name(feature)                       # Find gene name
            if name in genes:
                genes[name]["total"].append(rec.name)
                try:
                    if subfamily in genes[name]:
                        genes[name][subfamily].append(rec.name)
                    else:
                        genes[name][subfamily] = [rec.name]
                except:
                    print(f"{rec.name} subfamily undefined")
            else:
                genes[name] = {"total": [rec.name]}
                try:
                    genes[name][subfamily] = [rec.name]
                except:
                    print(f"{rec.name} subfamily undefined")
print(genes)

subfamilies = ['Agabinae', 'Colymbetinae', 'Copelatinae', 'Coptotominae', 'Cybistrinae', 'Dytiscinae', 'Hydrodytinae',
               'Hydroporinae', 'Laccophilinae', 'Lancetinae', 'Matinae']

# Write CSV metadata file
with open("metadata.csv", "w") as file:     # Open output file
    writer = csv.writer(file)               # Name writer object
    writer.writerow(
        ["Gene Name", "Count", 'Agabinae', 'Colymbetinae', 'Copelatinae', 'Coptotominae', 'Cybistrinae', 'Dytiscinae',
         'Hydrodytinae', 'Hydroporinae', 'Laccophilinae', 'Lancetinae', 'Matinae'])  # Write column names

for name, lists in genes.items():
    print(f"{len(lists['total'])} records for {name}")
    with open("metadata.csv", "a") as file:
        row = [name, len(lists["total"])]
        for sub in subfamilies:
            if sub in lists:
                row.append(len(lists[sub]))
            else:
                row.append('')
        writer = csv.writer(file)
        writer.writerow(row)
