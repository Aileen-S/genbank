#!/usr/bin/env python3

# python3 get_from_ids.py -e mixedupvoyage@gmail.com -f test.txt -i
# python3 get_from_ids.py -e mixedupvoyage@gmail.com -x -l txid2770180,txid2770181,txid2770182
# Get genbanks from list of GenBank ID numbers, accessions or taxon IDs.
# Have not added chuck search: only retrieves up to 10,000 records.


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
parser = argparse.ArgumentParser(description="Fetch metadata and fastas from specified GenBank accession/ID numbers. "
                                             "Input refs either in command with -r flag, or listed in text file using "
                                             "-f flag.", formatter_class=MultilineFormatter)
parser.add_argument("-l", "--list", type=str, help="GenBank ID/accession number(s). For multiple records, format is ref1,ref2,ref3.")
parser.add_argument("-f", "--file", type=str, help="Text file containing list of GenBank ID/accession refs, with one ref per line.")
parser.add_argument("-x", "--txid", action="store_true", help="Specify if input refs are NCBI taxon IDs.")
parser.add_argument('-i', '--fasta_id', action="store_true", help="Print taxon ID rather than accession in output fastas.")
parser.add_argument("-e", "--email", type=str, help="Your email registered with NCBI")

args = parser.parse_args()

#args = argparse.Namespace(file="test.txt", email='aileen.scott@nhm.ac.uk') # This is how I step through the script interactively
Entrez.email = args.email

# Gene name variants dict
genes = {"12S": ["12S", "12S RIBOSOMAL RNA", "12S RRNA", 'RRNS'],
         "16S": ["16S", "16S RIBOSOMAL RNA", "16S RRNA", "RRNL"],
         "ATP6": ['ATP SYNTHASE F0 SUBUNIT 6', 'APT6', 'ATP SYNTHASE A0 SUBUNIT 6', 'ATP SYNTHASE SUBUNIT 6', 'ATP SYNTHASE FO SUBUNIT 6', 'ATPASE6', 'ATPASE SUBUNIT 6', 'ATP6'],
         "ATP8": ['ATP SYNTHASE F0 SUBUNIT 8', 'APT8', 'ATP SYNTHASE A0 SUBUNIT 8', 'ATP SYNTHASE SUBUNIT 8', 'ATP SYNTHASE FO SUBUNIT 8', 'ATPASE8', 'ATPASE SUBUNIT 8', 'ATP8'],
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
         "18S": ["18S", "18S RIBOSOMAL RNA", "18S RRNA", "18S SMALL SUBUNIT RIBOSOMAL RNA", "SMALL SUBUNIT RIBOSOMAL RNA"],
         "28S": ["28S RIBOSOMAL RNA", "28S RRNA", "28S LARGE SUBUNIT RIBOSOMAL RNA", 'LARGE SUBUNIT RIBOSOMAL RNA'],
         "AK": ["AK", "ARGININE KINASE", "ARGK", "ARGKIN", "ARGS", "ARK"],
         "CAD": ["CAD", "CAD FRAGMENT 1", "CARBAMOYLPHOSPHATE SYNTHETASE"],
         "EF1A": ["EF1-ALPHA", "EF1A", "ELONGATION FACTOR 1 ALPHA", "ELONGATION FACTOR 1-ALPHA"],
         "H3": ["H3", "HISTONE 3", "HISTONE H3", "HIS3"],
         "RNApol": ["RNA POL II", "RNA POL2", "RNA POLYMERASE II LARGE SUBUNIT"],
         "Wg": ["WG", "WINGLESS", "WNG", "WNT", "WNT1", "WNT-4"]}


# Write CSV metadata file
with open("metadata.csv", "w") as file:     # Open output file
    writer = csv.writer(file)               # Name writer object
    writer.writerow(
        ["Accession", "Taxon ID", "Description", '18S', "28S", "AK", "CAD", 'EF1A', 'H3', 'RNApol', 'Wg',
         '12S', '16S', 'ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6',
         "Suborder", "Superfamily", "Family", "Subfamily", "Tribe", 'Genus', "Species", 'Label', "Date Late Modified",
         "Date Collected", "Country", "Region", "Lat/Long", "Ref1 Author", "Ref1 Title", "Ref1 Journal", "Ref2 Author",
         "Ref2 Title", "Ref2 Journal", "Ref3 Author", "Ref3 Title", "Ref3 Journal"])

gen = ['18S', '28S', 'AK', 'CAD', 'EF1A', 'H3', 'RNApol', 'Wg', '12S', '16S', 'ATP6', 'ATP8',
       'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']

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

if args.txid:
    # Get taxon IDs from command line input list format txid1,txid2,txid3
    if args.list:
        txids = args.list.split(",")
    # Get taxon IDs from file
    if args.file:
        if args.file:
            txids = []
            file = open(args.file)
            lines = file.readlines()
            for line in lines:
                txid = line.strip()
                txids.append(txid)

    # Search GenBank for txids and return list of accessions.
    ids = []
    for txid in txids:
        handle = Entrez.esearch(db="nucleotide", term=f"txid{txid}[Orgn]", rettype="gb", retmode="text",
                                retmax=10000)  # Get GenBanks
        record = Entrez.read(handle)
        ids = ids + record["IdList"]  # Get list of accessions
    id_str = ",".join(ids)

else:
    if args.list:
        id_str = args.list

    if args.file:
        ids = []
        file = open(args.file)
        lines = file.readlines()
        for line in lines:
            line.strip()
            ids.append(line)
        id_str = ",".join(ids)

# Fetch records from GenBank
handle = Entrez.efetch(db="nucleotide", id=id_str, rettype="gb", retmode="text")  # Get GenBanks
record = SeqIO.parse(handle, "gb")
sequences = []
for rec in record:
    x += 1
    spec = rec.annotations["organism"]
    specfasta = spec.replace(" ","_")
    taxonomy = rec.annotations["taxonomy"][10:15]
    taxonomy.extend([""] * (5 - len(taxonomy)))
    if taxonomy[4] == "Cybistrini":
        taxonomy[3] = "Cybistrinae"
    tax = f"_{taxonomy[2]}_{taxonomy[3]}_{taxonomy[4]}_{specfasta}"
    db_xref = rec.features[0].qualifiers["db_xref"]
    for ref in db_xref:
        if "taxon" in ref:                                  # Get NCBI taxon, rather than BOLD cross ref
            txid = "".join(filter(str.isdigit, ref))         # Extract numbers from NCBI taxon value
    gbid = rec.name
    if "country" in rec.features[0].qualifiers:
        location = rec.features[0].qualifiers["country"][0]
        if ":" in location:
            country, region = location.split(":")
        else:
            country = location
            region = ""
    else:
        country = ""
        region = ""
    if "lat_lon" in rec.features[0].qualifiers:
        latlon = rec.features[0].qualifiers["lat_lon"][0]
    else:
        latlon = ""
    if "collection_date" in rec.features[0].qualifiers:
        c_date = rec.features[0].qualifiers["collection_date"][0]
    else:
        c_date = ""
    refs = []
    for ref in rec.annotations["references"]:
        refs.append(ref.authors)
        refs.append(ref.title)
        refs.append(ref.journal)
    row = [gbid, txid, rec.description]          # Start row of metadata for CSV

    # Get sequences for fastas and csv
    feats = {}                                      # Gene length dict for metadata
    for feature in rec.features:
        type = feature.type
        if type not in ('CDS', 'rRNA', 'mRNA'):
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
        if stdname in sequencedict:                 # Save sequences for fasta
            if gbid not in sequencedict[stdname]:   # Avoid duplicate gene records
                sequencedict[stdname][gbid] = [txid, tax, seq]
        else:
            sequencedict[stdname] = {gbid: [txid, tax, seq]}
        if stdname not in feats.keys():             # Save lengths for metadata
            feats[stdname] = len(seq)


    # Continue row of metadata csv
    for g in gen:
        if g in feats:
            row.append(feats[g])
        else:
            row.append("")
    row.extend(taxonomy)
    gen_spec = spec.split(' ')
    genus = gen_spec[0]
    for k, v in subgenus.items():
        if genus in v:
            genus = k
    row2 = [genus, spec, gbid+tax, rec.annotations["date"], c_date, country, region, latlon,]
    row2.extend(refs)
    row = row + row2
    writer.writerow(row)


print(f"{str(x)} records found")
print(f"Unrecognised Genes {unrecgenes}")

for gene, records in sequencedict.items():
    file = open(f"{gene}.fasta", "w")
    for acc, rec in records.items():
        if args.fasta_id:
            f_id = rec[0]
        else:
            f_id = acc
        file.write(f">{f_id}{rec[1]}\n{rec[2]}\n")
    print(f"{len(records)} {gene} records saved to {gene}.fasta")
print("Metadata saved to metadata.csv")


