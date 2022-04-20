# Get list of species (db_xref) from GenBank and remove duplicates
# Each species is a dict
# For each species, each gene is a dict containing all the sequence lenghts
# Find longest, extract and save as fasta

import argparse
from Bio import Entrez
Entrez.email = "aileen.scott@nhm.ac.uk"


parser = argparse.ArgumentParser(description="Get list of species names of requested taxa.")
parser.add_argument("-t", "--taxon", type=str)  # Define command line inputs.
#parser.add_argument("-g", "--gene", type=str)   # -- means argument is optional, input with flags
#parser.add_argument("-l", "--length", type=str)
args = parser.parse_args()         # Process input args from command line

handle = Entrez.esearch(db="nucleotide", term=f"{args.taxon}", retmax=10)
record = Entrez.read(handle)
accs = record["IdList"]
print(str(record["Count"]) + " records found")

accs_str = ",".join(accs)
handle = Entrez.esummary(db="nucleotide", id=accs_str)
record = Entrez.read(handle)
TaxList = []
Accessions = []
for rec in record:
    #print(rec["TaxId"])
    TaxId = rec["TaxId"]        # Look at db_xref taxon identifiers
    if TaxId in TaxList:        # Check if taxon is already in the list
        continue                # If so, skip
    else:
        TaxList.append(TaxId)           # If not, add to list
        Accessions.append(rec["Id"])    # Save accession
print(str(len(TaxList)) + " unique species saved")


