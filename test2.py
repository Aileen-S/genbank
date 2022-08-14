# python3 get_from_ids.py -e mixedupvoyage@gmail.com -f test.txt -x
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


Entrez.email = 'mixedupvoyage@gmail.com'

# Gene name variants dict
genes = {"12S": ["12S", "12S RIBOSOMAL RNA", "12S RRNA"],
         "16S": ["16S", "16S RIBOSOMAL RNA", "16S RRNA"],
         "18S": ["18S", "18S RIBOSOMAL RNA", "18S RRNA", "18S SMALL SUBUNIT RIBOSOMAL RNA"],
         "EF1A": ["EF1-ALPHA", "EF1A", "ELONGATION FACTOR 1 ALPHA", "ELONGATION FACTOR 1-ALPHA", "EF-1A"],
         "H3": ["H3", "HISTONE 3", "HISTONE H3", "HIS3"],
         "Wg": ["WG", "WINGLESS", "WNG", "WNT", "WNT1", "WNT-4"],
         "ATP6": ['ATP SYNTHASE F0 SUBUNIT 6', 'APT6', 'ATP SYNTHASE A0 SUBUNIT 6', 'ATP SYNTHASE SUBUNIT 6', 'ATP SYNTHASE FO SUBUNIT 6', 'ATPASE6', 'ATPASE SUBUNIT 6', 'ATP6'],
         "ATP8": ['ATP SYNTHASE F0 SUBUNIT 8', 'APT8', 'ATP SYNTHASE A0 SUBUNIT 8', 'ATP SYNTHASE SUBUNIT 8', 'ATP SYNTHASE FO SUBUNIT 8', 'ATPASE8', 'ATPASE SUBUNIT 8', 'ATP8'],
         "COX1": ['CYTOCHROME C OXIDASE SUBUNIT 1', 'CYTOCHROME OXIDASE SUBUNIT I', 'CYTOCHROME C OXIDASE SUBUNIT I', 'COXI', 'CO1', 'COI', 'CYTOCHROME COXIDASE SUBUNIT I', 'CYTOCHROME OXIDASE SUBUNIT 1', 'CYTOCHROME OXYDASE SUBUNIT 1', 'COX1'],
         "COX2": ['CYTOCHROME C OXIDASE SUBUNIT 2', 'CYTOCHROME OXIDASE SUBUNIT II', 'CYTOCHROME C OXIDASE SUBUNIT II', 'COXII', 'CO2', 'COII', 'CYTOCHROME COXIDASE SUBUNIT II', 'CYTOCHROME OXIDASE SUBUNIT 2', 'CYTOCHROME OXYDASE SUBUNIT 2', 'COX2'],
         "COX3": ['CYTOCHROME C OXIDASE SUBUNIT 3', 'CYTOCHROME OXIDASE SUBUNIT III', 'CYTOCHROME C OXIDASE SUBUNIT III', 'COXII', 'CO3', 'COIII', 'CYTOCHROME COXIDASE SUBUNIT III', 'CYTOCHROME OXIDASE SUBUNIT 3', 'CYTOCHROME OXYDASE SUBUNIT 3', 'COX3'],
         "CYTB": ['CYTOCHROME B', 'CYB', 'COB', 'COB / CYTB', 'CYTB'],
         "ND1": ['NAD1', 'NSD1', 'NADH1', 'NADH DEHYDROGENASE SUBUNIT I', 'NADH DEHYDROGENASE SUBUNIT 1', 'NADH DESHYDROGENASE SUBUNIT 1', 'NAD1-0', 'ND1'],
         "ND2": ['NAD2', 'NSD2', 'NADH2', 'NADH DEHYDROGENASE SUBUNIT II', 'NADH DEHYDROGENASE SUBUNIT 2', 'NADH DESHYDROGENASE SUBUNIT 2', 'NAD2-0', 'ND2'],
         "ND3": ['NAD3', 'NSD3', 'NADH3', 'NADH DEHYDROGENASE SUBUNIT III', 'NADH DEHYDROGENASE SUBUNIT 3', 'NADH DESHYDROGENASE SUBUNIT 3', 'NAD3-0', 'ND3'],
         "ND4": ['NAD4', 'NSD4', 'NADH4', 'NADH DEHYDROGENASE SUBUNIT IV', 'NADH DEHYDROGENASE SUBUNIT 4', 'NADH DESHYDROGENASE SUBUNIT 4', 'NAD4-0', 'ND4'],
         "ND4L": ['NAD4L', 'NSD4L', 'NADH4L', 'NADH DEHYDROGENASE SUBUNIT IVL', 'NADH DEHYDROGENASE SUBUNIT 4L', 'NADH DESHYDROGENASE SUBUNIT 4L', 'NAD4L-0', 'ND4L'],
         "ND5": ['NAD5', 'NSD5', 'NADH5', 'NADH DEHYDROGENASE SUBUNIT V', 'NADH DEHYDROGENASE SUBUNIT 5', 'NADH DESHYDROGENASE SUBUNIT 5', 'NAD5-0', 'ND5'],
         "ND6": ['NAD6', 'NSD6', 'NADH6', 'NADH DEHYDROGENASE SUBUNIT VI', 'NADH DEHYDROGENASE SUBUNIT 6', 'NADH DESHYDROGENASE SUBUNIT 6', 'NAD6-0', 'ND6'],
         "28S": ["28S RIBOSOMAL RNA", "28S RRNA", "28S LARGE SUBUNIT RIBOSOMAL RNA"],
         "AK": ["AK", "ARGININE KINASE", "ARGK", "ARGKIN", "ARGS", "ARK"],
         "CAD": ["CAD", "CAD FRAGMENT 1"]}

gen = ['12S', '16S', '18S', 'EF1A', 'H3', 'Wg', 'ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']


# Write CSV metadata file
with open("metadata.csv", "w") as file:     # Open output file
    writer = csv.writer(file)               # Name writer object
    writer.writerow(
        ["Accession", "Taxon ID", "Description", "Gene Count",
         '12S Length', '16S Length', '18S Length', 'EF1A Length', 'H3 Length', 'Wg Length', 'ATP6 Length',
         'ATP8 Length', 'COX1 Length', 'COX2 Length', 'COX3 Length', 'CYTB Length', 'ND1 Length', 'ND2 Length',
         'ND3 Length', 'ND4 Length', 'ND4L Length', 'ND5 Length', 'ND6 Length',
         "Domain", "Kingdom", "Superphylum", "Phylum",
         "Subphylum", "Class", "Subclass", "Infraclass", "Superorder", "Order", "Suborder", "Superfamily", "Family",
         "Subfamily", "Tribe", "Species", "Country", "Region", "Lat/Long", "Date Late Modified", "Date Collected",
         "Ref1 Author", "Ref1 Title", "Ref1 Journal", "Ref2 Author", "Ref2 Title", "Ref2 Journal", "Ref3 Author",
         "Ref3 Title", "Ref3 Journal"])

unrecgenes = set()
sequencedict = {}
file = open("metadata.csv", "a")
writer = csv.writer(file)
x = 0  # Count records added to species dict.



file = open("sub_concat.gb")
record = SeqIO.parse(file, "gb")
sequences = []
for rec in record:
    x += 1
    spec = rec.annotations["organism"]
    taxonomy = rec.annotations["taxonomy"][0:15]
    db_xref = rec.features[0].qualifiers["db_xref"]
    for ref in db_xref:
        if "taxon" in ref:                                  # Get NCBI taxon, rather than BOLD cross ref
            tax = "".join(filter(str.isdigit, ref))         # Extract numbers from NCBI taxon value
    rec_id = rec.annotations["taxonomy"][13]
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
    row1 = [rec_id, tax, rec.description]          # Start row of metadata for CSV

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
            if rec_id not in sequencedict[stdname]:   # Avoid duplicate gene records
                sequencedict[stdname][rec_id] = seq
        else:
            sequencedict[stdname] = {rec_id: seq}
        if stdname not in feats.keys():             # Save lengths for metadata
            feats[stdname] = len(seq)


    # Continue row of metadata csv
    row1.append(len(feats))
    for g in gen:
        if g in feats:
            row1.append(feats[g])
        else:
            row1.append("")
    taxonomy.extend([""] * (15 - len(taxonomy)))
    for t in taxonomy:
        row1.append(t)
    row2 = [spec, country, region, latlon, rec.annotations["date"], c_date]
    for r in refs:
        row2.append(r)
    row = row1 + row2
    writer.writerow(row)


print(f"{str(x)} records found")
print(f"Unrecognised Genes {unrecgenes}")

for gene, ids in sequencedict.items():
    file = open(f"{gene}.fasta", "w")
    for rec_id, seq in ids.items():
        file.write(f">{rec_id}\n{seq}\n")
    print(f"{len(ids)} {gene} records saved to {gene}.fasta")
print("Metadata saved to metadata.csv")


