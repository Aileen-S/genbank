# Trying to use argparse. Input gene, taxon and desired minimum length.
# Search first for taxon, then parse XML file for gene using Thomas's genenamevarients
# Find desired gene
# Filter minimum gene length
# Choose one for each species
# Download genbank
# Extract sequences
import Bio.Entrez
from Bio import SeqIO
from Bio import Entrez
Entrez.email = "aileen.scott@nhm.ac.uk"
import argparse
import urllib

parser = argparse.ArgumentParser(description="Extract genes of specified taxon and minimum length from GenBank seach. Output as fasta.")
parser.add_argument("-t", "--taxon", type=str)  # Define command line inputs.
parser.add_argument("-g", "--gene", type=str)   # -- means argument is optional, input with flags
parser.add_argument("-l", "--length", type=str)
args = parser.parse_args()         # Process input args from command line

# GenBank search for taxon and gene
handle = Entrez.esearch(db="nuccore", term=f"{args.taxon}", retmax=10) # Make term argument into single string
record = Entrez.read(handle)
print(str(record["Count"]) + "record found")

# Use accession numbers to download GenBank records
# Couldn't work out how to use retstart to return >10,000 (for COX1)

# Get XML file
accessions = record["IdList"]
accessions_str = ",".join(accessions)
handle = Entrez.efetch(db="nucleotide", id=accessions_str, retmode="xml")
record = Bio.Entrez.read(handle)
#print(record, file=open("test.xml", "w"))
for rec in record:
    results = rec['GBSeq_feature-table'][1]['GBFeature_quals']
    print(results)
# Not sure what to call for to get gene names.
# Not sure there is enough info in xml file to get the data I need anyway, might need full GenBank file.

#print(record, file=open("test.xml", "w"))



