# AIMS
# Get list of species (db_xref) from GenBank and remove duplicates
# Make dict: species = {taxonID {gene {records}}}
# Find longest, extract and save as fasta.
# Separate fasta file for each gene, ready for alignment.

# PROBLEMS
# Lines from "for rec in record:" : unsure if loop syntax is right. Not saving all sequences.
# Better way to get taxonomy for CSV?
# Outdated records. Eg Dytiscinae split into Dytiscinae and Cybistrinae now. Need to manually check each?
# Possible to add filter to prioritise records with multiple genes, when filtering for length? Better for alignment.
# Alignment very gappy. Possible outliers, need to find these.
# Duplicate records in CSV file

# TO DO
# Need to add argparse option to search for specific genes
# Set min sequence length?
# Set max sequence length 20000[SLEN]
# How to include ATP8
# How to include COX1: USEARCH or manual split after alignment? Filter full length sequences to align separately?
# Subspecies problem
# Add nuclear genes/16S and name variants


#python3 get_species_list.py -e aileen.scott@nhm.ac.uk -t Agabus -m


import argparse
import urllib
import csv
from Bio import Entrez
from Bio import SeqIO
import textwrap as _textwrap
from collections import defaultdict

# Function definitions

# Opens url,
def loadnamevariants():
    conversion = {}
    url = "https://raw.githubusercontent.com/tjcreedy/constants/master/gene_name_variants.txt"
    fullparse = {}
    alltypes = set()
    for line in urllib.request.urlopen(url):
        line = line.decode('utf-8').strip()
        description, variants = line.split(":")
        name, annotype, fullname = description.split(";")
        variants = variants.split(',')
        variants.extend([name, fullname.upper()])

        fullvariants = []
        for v in [name] + variants:
            for g in ['', ' ']:
                v = v.replace(g, '')
                for s in ['', ' GENE', ' '+annotype.upper()]:
                    fullvariants.append(v+s)
                    conversion[v+s] = name

        alltypes.add(annotype)
        fullparse[name] = {'type': annotype, 'variants': fullvariants, 'product': fullname}
    return conversion, alltypes, fullparse


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


# Filter sequence length
# Returns true if sequence is in target range
def seqlength(x):
    url = "https://raw.githubusercontent.com/Aileen-S/genbank/main/sequencelength.txt?token=GHSAT0AAAAAABQZX632I3W6OLMRN2TREMRQYT3X2MQ"
    filter = {}
    for line in urllib.request.urlopen(url):
        line = line.decode('utf-8').strip()
        gene, lengths = line.split(":")
        actual, range = lengths.split(';')
        min, max = range.split(",")
        filter[gene] = [min, max]
    for k, v in filter.items():
        if str(x) in k:
            if int(v[0]) < x < int(v[1]):
                return True
            else:
                return False

# Write CSV metadata file
with open("metadata.csv", "w") as file:     # Open output file
    writer = csv.writer(file)               # Name writer object
    writer.writerow(
        ["Accession", "Taxon ID", "Species", "Gene", "Sequence Length", "Domain", "Kingdom", "Superphylum", "Phylum", "Subphylum", "Class", "Subclass",
         "Infraclass", "Superorder", "Order", "Suborder", "Superfamily", "Family", "Subfamily", "Tribe", "1", "2",
         "3"])  # Write column names


# Write row of metadata file
def writecsv(x):                                # x = genbank record
    gbid  = x.name                              # Get accessions (record.name gives accession, record.id gives version)
    spe = x.annotations["organism"]             # Get genus/species
    taxonomy = x.annotations["taxonomy"]        # Get higher taxonomy list
    row = [gbid, tax, spe, stdname, len(sequence)]                           # New list with accession and species names
    for level in taxonomy:
        row.append(level)                       # Add taxon levels to list
    with open("metadata.csv", "a") as file:
        writer = csv.writer(file)               # Name writer object
        writer.writerow(row)                    # Write row


# This reclasses the argparse.HelpFormatter object to have newlines in the help text for paragraphs
class MultilineFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(' ', text).strip()
        paragraphs = text.split('|n ')
        multiline_text = ''
        for paragraph in paragraphs:
            formatted_paragraph = _textwrap.fill(paragraph, width, initial_indent=indent,
                                                 subsequent_indent=indent
                                                 ) + '\n\n'
            multiline_text = multiline_text + formatted_paragraph
        return multiline_text


# Argument parser
parser = argparse.ArgumentParser(description="Search GenBank, retrive gene sequences and save as fasta.", formatter_class=MultilineFormatter)
parser.add_argument("-t", "--taxon", type=str, help="Taxon of interest")
parser.add_argument("-m", "--mpc", action="store_true", help="Save only mitochondrial protein coding genes")
#parser.add_argument("-g", "--gene", type=str, help="Gene(s) of interest: if not specified, all genes will be retrieved")
parser.add_argument("-e", "--email", type=str, help="Your email registered with NCBI")
parser.add_argument("-n", "--nuclear", action="store_true", help="Search for nuclear as well as mitochondrial genes.")
# parser.add_argument("-l", "--length", type=str)


# Start the actual script

args = parser.parse_args()         # Process input args from command line
# args = parser.parse_args('-t Dytiscidae -e '.split(' ')) # This is how I step through the script interactively

# Get name variants
nameconvert, types, namevariants = loadnamevariants()


# Set up for unrecognised genes
#unrecgenes = defaultdict(list)
unrecgenes = set()

# Set for mitochondiral non-coding genes (filtered out if -m argument used)
mncgenes = set()

Entrez.email = args.email

# Generate search term to get all sequences in the search taxonomy
# - if -n option not used, then include "mitochondrial" in search term.
basesearch = f"(\"{args.taxon}\"[Organism] OR \"{args.taxon}\"[All Fields])"\
             f"{'' if args.nuclear else ' AND mitochondrion[filter]'}"

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

mpc = ["ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"]
x = 0
y = 0
species = {}
for tax in taxids:
    y += 1
    if y % 100 == 0:
        print(f"Downloading GenBank records for taxon IDs {y+1} to {y+100}" if (y+100) > len(taxids) else
              f"Downloading GenBank records for taxon IDs {y+1} to {len(taxids)}")
    handle = Entrez.esearch(db="nucleotide", term=f"txid{tax}")       # Search for all records for each taxon id
    record = Entrez.read(handle)
    accs   = record["IdList"]                                         # Get accessions
    accstr = ",".join(accs)                                           # Join into string for efetch
    handle = Entrez.efetch(db="nucleotide", id=accstr, rettype="gb", retmode="text")  # Get GenBanks
    record = SeqIO.parse(handle, "gb")
    for rec in record:
        for feature in rec.features:
            type = feature.type
            if type not in ('CDS', 'rRNA'):
                continue  # skip the rest of the current iteration of this loop
            name = get_feat_name(feature)                       # Use function to search for gene names
            if name in nameconvert:
                stdname = nameconvert[name]                    # If gene name in namevariants, convert to standard name
                if args.mpc:
                    if stdname in mpc:   # Filter: keep only mitochondrial protein coding genes if -m argument used
                        pass
                    else:
                        mncgenes.add(stdname)
                        continue
                sequence = rec[feature.location.start:feature.location.end]
                output = [stdname, rec.name, rec.description, type, len(sequence), str(sequence.seq)]
                if seqlength(len(sequence)) == True:
                    if tax in species:                              # If taxon ID in dict
                        if stdname in species[tax]:                 # If gene in dict for that taxon ID
                            species[tax][stdname].append(output)    # Add gene info list to dict
                            x += 1
                        else:
                            species[tax][stdname] = [output]      # Otherwise add to dict with new key
                        x += 1
                    else:
                        species[tax] = {stdname: [output]}      # Otherwise add to dict with new key
                        x += 1
                else:
                    unrecgenes.add(name)

print(f"{str(x)} gene records saved to species dict")
#print(species)

if args.mpc:
    print("\nMitochondrial non-coding genes, not saved:")
    print(mncgenes)

print("\nUnrecognised Genes")
print(unrecgenes)

# Temporary: write dict to file to play with length filter
#import json
#with open("testspeciesdict.txt", "w") as file:
    #file.write(json.dumps(species))


# Set first record as max value, iterate through records and replace if another sequence is longer.
# Record format = [gene, ID, rec.description, type, length, sequence]
def findmax(x):
    maxrec = x[0]
    for record in x:
        if record[4] > maxrec[4]:
            maxrec = record
    return maxrec


# Dict for longest sequences, key is gene stdname, value is list of records
longest = {}
for tax, stdname in species.items():
    for gene, list in stdname.items():
        new = findmax(list)
        if gene in longest:
            longest[gene].append(new)
        else:
            longest[gene] = [new]
#print(species)
#print(longest)

# Save each gene list to separate fasta file
for gene, records in longest.items():
    file = open(f"{gene}test.fasta", "w")
    for rec in records:
        writecsv(rec)
        file.write(">" + rec[1] + " " + rec[2] + "\n" + rec[5] + "\n")



