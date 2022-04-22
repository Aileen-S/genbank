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

handle = Entrez.esearch(db="nucleotide", term=f"{args.taxon}", retmax=3)# Search for all records of specified taxon
record = Entrez.read(handle)
accs   = record["IdList"]                                               # Save accession numbers
print(str(record["Count"]) + " records found")                          # Print total records
print(record)

accs_str = ",".join(accs)
handle = Entrez.esummary(db="nucleotide", id=accs_str)                  # Get esummary for accessions
record = Entrez.read(handle)
TaxSet = set()
for rec in record:                                                      # Save taxon ids in a set (removes duplicates)
    Tax = ("txid"+str(int(rec["TaxId"])))                               # Add "txid before number for esearch
    TaxSet.add(Tax)

print(str(len(TaxSet)) + " unique species saved")                       # Print total taxon ids

Species = {}
for Tax in TaxSet:
    handle = Entrez.esearch(db="nucleotide", term=Tax)                  # Search for all records for each taxon id
    record = Entrez.read(handle)
    IdList = record["IdList"]                                           # Get accessions
    IdStr  = ",".join(IdList)                                           # Join into string for efetch
    handle = Entrez.efetch(db="nucleotide", id=IdStr, rettype="gb")     # Get GenBank record for accession list
    record = handle.read()

# Need to make dict with species, genes and sequences
# Work out name variants

#for k,v in Species.items():

# Then efetch from GIs
# Add mitochondrion to term
# Max sequence length 20000[SLEN]