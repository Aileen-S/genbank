#!/usr/bin/env python3


# python3 ptp_get_metadata.py -e mixedupvoyage@gmail.com -f test.txt

import argparse
import csv
from Bio import Entrez
from Bio import SeqIO
import textwrap as _textwrap


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
parser = argparse.ArgumentParser(description="Combine metadata from multiple genes")
parser.add_argument("-f", "--file", type=str, help="Text file containing list of GenBank ID/accession refs, with one ref per line.")
parser.add_argument("-e", "--email", type=str, help="Your email registered with NCBI")

args = parser.parse_args()

#args = argparse.Namespace(file="test.txt", email='aileen.scott@nhm.ac.uk') # This is how I step through the script interactively
Entrez.email = args.email

# Gene name variants dict
genes = {"12S": ["12S", "12S RIBOSOMAL RNA", "12S RRNA", 'RRNS'],
         "16S": ["16S", "16S RIBOSOMAL RNA", "16S RRNA", "RRNL"],
         "18S": ["18S", "18S RIBOSOMAL RNA", "18S RRNA", "18S SMALL SUBUNIT RIBOSOMAL RNA", "SMALL SUBUNIT RIBOSOMAL RNA"],
         "ATP6": ['ATP SYNTHASE F0 SUBUNIT 6', 'APT6', 'ATP SYNTHASE A0 SUBUNIT 6', 'ATP SYNTHASE SUBUNIT 6', 'ATP SYNTHASE FO SUBUNIT 6', 'ATPASE6', 'ATPASE SUBUNIT 6', 'ATP6'],
         "COX1": ['CYTOCHROME C OXIDASE SUBUNIT 1', 'CYTOCHROME OXIDASE SUBUNIT I', 'CYTOCHROME C OXIDASE SUBUNIT I', 'COXI', 'CO1', 'COI', 'CYTOCHROME COXIDASE SUBUNIT I', 'CYTOCHROME OXIDASE SUBUNIT 1', 'CYTOCHROME OXYDASE SUBUNIT 1', 'COX1'],
         "COX2": ['CYTOCHROME C OXIDASE SUBUNIT 2', 'CYTOCHROME OXIDASE SUBUNIT II', 'CYTOCHROME C OXIDASE SUBUNIT II', 'COXII', 'CO2', 'COII', 'CYTOCHROME COXIDASE SUBUNIT II', 'CYTOCHROME OXIDASE SUBUNIT 2', 'CYTOCHROME OXYDASE SUBUNIT 2', 'COX2'],
         "COX3": ['CYTOCHROME C OXIDASE SUBUNIT 3', 'CYTOCHROME OXIDASE SUBUNIT III', 'CYTOCHROME C OXIDASE SUBUNIT III', 'COXII', 'CO3', 'COIII', 'CYTOCHROME COXIDASE SUBUNIT III', 'CYTOCHROME OXIDASE SUBUNIT 3', 'CYTOCHROME OXYDASE SUBUNIT 3', 'COX3'],
         "CYTB": ['CYTOCHROME B', 'CYB', 'COB', 'COB / CYTB', 'CYTB', "COB/CYTB"],
         "ND1": ['NAD1', 'NSD1', 'NADH1', 'NADH DEHYDROGENASE SUBUNIT I', 'NADH DEHYDROGENASE SUBUNIT 1', 'NADH DESHYDROGENASE SUBUNIT 1', 'NAD1-0', 'ND1'],
         "ND2": ['NAD2', 'NSD2', 'NADH2', 'NADH DEHYDROGENASE SUBUNIT II', 'NADH DEHYDROGENASE SUBUNIT 2', 'NADH DESHYDROGENASE SUBUNIT 2', 'NAD2-0', 'ND2'],
         "ND3": ['NAD3', 'NSD3', 'NADH3', 'NADH DEHYDROGENASE SUBUNIT III', 'NADH DEHYDROGENASE SUBUNIT 3', 'NADH DESHYDROGENASE SUBUNIT 3', 'NAD3-0', 'ND3'],
         "ND4": ['NAD4', 'NSD4', 'NADH4', 'NADH DEHYDROGENASE SUBUNIT IV', 'NADH DEHYDROGENASE SUBUNIT 4', 'NADH DESHYDROGENASE SUBUNIT 4', 'NAD4-0', 'ND4'],
         "ND4L": ['NAD4L', 'NSD4L', 'NADH4L', 'NADH DEHYDROGENASE SUBUNIT IVL', 'NADH DEHYDROGENASE SUBUNIT 4L', 'NADH DESHYDROGENASE SUBUNIT 4L', 'NAD4L-0', 'ND4L'],
         "ND5": ['NAD5', 'NSD5', 'NADH5', 'NADH DEHYDROGENASE SUBUNIT V', 'NADH DEHYDROGENASE SUBUNIT 5', 'NADH DESHYDROGENASE SUBUNIT 5', 'NAD5-0', 'ND5'],
         "ND6": ['NAD6', 'NSD6', 'NADH6', 'NADH DEHYDROGENASE SUBUNIT VI', 'NADH DEHYDROGENASE SUBUNIT 6', 'NADH DESHYDROGENASE SUBUNIT 6', 'NAD6-0', 'ND6'],
         "H3": ["H3", "HISTONE 3", "HISTONE H3", "HIS3"],
         "Wg": ["WG", "WINGLESS", "WNG", "WNT", "WNT1", "WNT-4"]}

#"ATP8": ['ATP SYNTHASE F0 SUBUNIT 8', 'APT8', 'ATP SYNTHASE A0 SUBUNIT 8', 'ATP SYNTHASE SUBUNIT 8', 'ATP SYNTHASE FO SUBUNIT 8', 'ATPASE8', 'ATPASE SUBUNIT 8', 'ATP8'],
#"RNApol": ["RNA POL II", "RNA POL2", "RNA POLYMERASE II LARGE SUBUNIT"],
#"28S": ["28S RIBOSOMAL RNA", "28S RRNA", "28S LARGE SUBUNIT RIBOSOMAL RNA", 'LARGE SUBUNIT RIBOSOMAL RNA'],
#"AK": ["AK", "ARGININE KINASE", "ARGK", "ARGKIN", "ARGS", "ARK"],
#"CAD": ["CAD", "CAD FRAGMENT 1", "CARBAMOYLPHOSPHATE SYNTHETASE"],
#"EF1A": ["EF1-ALPHA", "EF1A", "ELONGATION FACTOR 1 ALPHA", "ELONGATION FACTOR 1-ALPHA"],

# Write CSV metadata file
with open("metadata.csv", "w") as file:     # Open output file
    writer = csv.writer(file)               # Name writer object
    writer.writerow(
        ["Accessions", "Taxon ID", "Species", '12S', "16S", "18S", "H3", 'Wg',
         'ATP6', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6',
         "Suborder", "Superfamily", "Family", "Subfamily", "Tribe", 'Genus'])

gen = ['12S', '16S', '18S', 'H3', 'Wg', 'ATP6',
       'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']

cds = ['ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6', 'AK', 'CAD', 'EF1A', 'H3', 'RNApol', 'Wg']


subgenus = {'Agabus': ['Acatodes', 'Gaurodytes'],
            'Platynectes': ['Agametrus', 'Australonectes', 'Gueorguievtes', 'Leuronectes'],
            'Cybister': ['Megadytoides', 'Melanectes', 'Neocybister'],
            'Megadytes': ['Bifurcitus', 'Paramegadytes', 'Trifurcitus'],
            'Acilus': ['Homoeolytrus'],
            'Hydaticus': ['Prodaticus'],
            'Clypeodytes': ['Hypoclypeus', 'Paraclypeus'],
            'Clemnius': ['Cyclopius'],
            'Hygrotus': ['Coelambus', 'Heroceras', 'Herophydrus', 'Hyphoporus', 'Leptolambus'],
            'Rhantus': ['Anisomera', 'Senilites'],
            'Aglymbus': ['Rugosus'],
            'Exocelina': ['Papuadytes'],
            'Paroster': ['Terradessus']}

unrecgenes = set()
sequencedict = {}
file = open("metadata.csv", "a")
writer = csv.writer(file)
x = 0  # Count records added to species dict.

# Get IDs from argparse input

ids = []
file = open(args.file)
lines = file.readlines()
for line in lines:
    line.strip()
    ids.append(line)
id_str = ",".join(ids)


txids = {}
# Fetch records from GenBank
handle = Entrez.efetch(db="nucleotide", id=id_str, rettype="gb", retmode="text")  # Get GenBanks
record = SeqIO.parse(handle, "gb")
sequences = []
for rec in record:
    x += 1
    y = 0
    species = rec.annotations["organism"]
    genus_spec = species.split(' ', 1)   #.split('mango', 1)[1]
    for k, v in subgenus.items():
        if genus_spec[0] in v:
            genus_spec[0] = k
    genus_spec[1] = genus_spec[1].replace(" ", "_")
    taxonomy = rec.annotations["taxonomy"][10:15]
    taxonomy.extend([""] * (5 - len(taxonomy)))
    if taxonomy[4] == "Cybistrini":
        taxonomy[3] = "Cybistrinae"
    tax = f"_{taxonomy[2]}_{taxonomy[3]}_{taxonomy[4]}_{genus_spec[0]}_{genus_spec[1]}"
    db_xref = rec.features[0].qualifiers["db_xref"]
    for ref in db_xref:
        if "taxon" in ref:                                  # Get NCBI taxon, rather than BOLD cross ref
            txid = "".join(filter(str.isdigit, ref))         # Extract numbers from NCBI taxon value
    gbid = rec.name
    if txid in txids.keys():
        txids[txid]['gbids'].append(gbid)
    else:
        txids[txid] = {'gbids': [gbid],
                       'species': species,
                       'taxonomy': taxonomy,
                       'gs': genus_spec}


    # Get sequences for fastas and csv
    feats = {}                                      # Gene length dict for metadata
    for feature in rec.features:
        type = feature.type
        if type not in ('CDS', 'rRNA'):
            continue  # skip to next feature
        name = get_feat_name(feature)
        stdname = ""
        for k, v in genes.items():
            if name in v:
                stdname = k
        if stdname == "":
            unrecgenes.add(name)
            continue
        else:
            seq = feature.extract(rec.seq)
            y += 1
            txids[txid][stdname] = len(seq)

for txid, data in txids.items():
    row = [",".join(data['gbids']), txid, data['species']]          # Start row of metadata for CSV

# Continue row of metadata csv
    for g in gen:
        if g in data.keys():
            row.append(data[g])
        else:
            row.append("")
    row.extend(data['taxonomy'])
    row.append(data['gs'][0])
    writer.writerow(row)

print(f"{str(x)} records found")
print("Metadata saved to metadata.csv")
