# Get list of species (db_xref) from GenBank and remove duplicates
# Each species is a dict
# For each species, each gene is a dict containing all the sequence lenghts
# Find longest, extract and save as fasta

import argparse
import urllib
from Bio import Entrez
from Bio import SeqIO
from collections import defaultdict
Entrez.email = "aileen.scott@nhm.ac.uk"


def loadnamevariants():
    output = {}
    url = "https://raw.githubusercontent.com/tjcreedy/genenames/master/gene_name_variants.txt"
    for line in urllib.request.urlopen(url):
        line = line.decode('utf-8').strip()
        name = line.split(";")[0]
        annotype = line.split(":")[0].split(";")[1]
        variants = line.split(":")[1].split(",")
        for v in variants:
            for g in ['', ' ']:
                v = v.replace(g, '')
                for s in ['',' GENE', ' '+annotype.upper()]:
                    output[v+s] = name
    return(output)


namevariants = loadnamevariants()

# Set up for unrecognised genes
unrecgenes = defaultdict(list)


# Get gene name from record.features
def get_feat_name(feat):
    featname = "unknown"
    nametags = ['gene', 'product', 'label', 'standard_name']  # Search these four keys for gene name
    if any(t in feat.qualifiers.keys() for t in nametags):
        for t in nametags:
            if t in feat.qualifiers.keys():
                featname = feat.qualifiers[t][0].upper()
                break
    return featname


# Change gene name to standard name
def set_feat_name(feat, name):
    nametags = ['gene', 'product', 'label', 'standard_name']
    if any(t in feat.qualifiers.keys() for t in nametags):
        for t in nametags:
            if t in feat.qualifiers.keys() :
                feat.qualifiers[t][0] = name
    return feat


parser = argparse.ArgumentParser(description="Get list of species names of requested taxa.")
parser.add_argument("-t", "--taxon", type=str)  # Define command line inputs.
parser.add_argument("-g", "--gene", type=str)   # -- means argument is optional, input with flags
# parser.add_argument("-l", "--length", type=str)
args = parser.parse_args()         # Process input args from command line

handle = Entrez.esearch(db="nucleotide", term=f"{args.taxon}", retmax=10)  # Search for all records of specified taxon
record = Entrez.read(handle)
accs   = record["IdList"]                                               # Save accession numbers
print(str(record["Count"]) + " records found")                          # Print total records

accstr = ",".join(accs)
handle = Entrez.esummary(db="nucleotide", id=accstr)                  # Get esummary for accessions
record = Entrez.read(handle)
taxids = set()                                                          # Make empty set to save taxon ids
for rec in record:                                                      # Save taxon ids in a set (removes duplicates)
    tax = ("txid"+str(int(rec["TaxId"])))                               # Add "txid" before id number for esearch
    taxids.add(tax)

print(str(len(taxids)) + " unique species saved")                       # Print total taxon ids

species = {}
for tax in taxids:
    handle = Entrez.esearch(db="nucleotide", term=tax)                # Search for all records for each taxon id
    record = Entrez.read(handle)
    accs   = record["IdList"]                                         # Get accessions
    accstr = ",".join(accs)                                           # Join into string for efetch
    handle = Entrez.efetch(db="nucleotide", id=accstr, rettype="gb", retmode="text")  # Get GenBanks
    record = SeqIO.parse(handle, "gb")
    for rec in record:
        for feature in rec.features:
            type = feature.type                                 # Retrieve the feature type, might be useful later
            name = get_feat_name(feature)                       # Use function to search for gene names
            if name in namevariants:
                stdname = namevariants[name]                    # If gene name in namevariants, convert to standard name
                sequence = rec[feature.location.start:feature.location.end]
                output = [rec.name, type, len(sequence)]        # Save relevent info in list (Add sequence when working)
                if tax in species:                              # If taxon ID in dict
                    if stdname in species[tax]:                 # If gene in dict for that taxon ID
                        species[tax][stdname].append(output)    # Add gene info list to dict
                    else:
                        species[tax] = {stdname: [output]}          # Add stdname as new key
                else:
                    species[tax] = {stdname: [output]}      # Add tax as new key

            else:
                unrecgenes[name].append(name)            # If gene name not in namevarants, save to list to check later

print("\nSpecies Dict")
print(species)
print("\nUnrecognised Genes")
print(unrecgenes)




# Fix gene name output
# Work out name variants

# for k,v in Species.items():

# Add mitochondrion to term
# Max sequence length 20000[SLEN]
