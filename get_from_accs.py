

#python3 get_other_genes.py -e mixedupvoyage@gmail.com -t Eretes


import argparse
import csv
from Bio import Entrez
from Bio import SeqIO
import textwrap as _textwrap
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
Entrez.email = "mixedupvoyage@gmail.com"

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

# Get record/feature dict from genbank record
def get_output(rec):
    sequence = rec[feature.location.start:feature.location.end]
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
    output = {"gene": stdname,
              "gbid": rec.name,
              "txid": tax,
              "description": rec.description,
              "spec": rec.annotations["organism"],
              "rec date": rec.annotations["date"],
              "c date": c_date,
              "taxonomy": rec.annotations["taxonomy"][0:15],
              "type": type,
              "length": len(sequence),
              "seq": feature.extract(rec.seq),
              "country": country,
              "region": region,
              "latlon": latlon,
              "refs": refs}
    return output


# Write CSV metadata file
with open("metadata.csv", "w") as file:     # Open output file
    writer = csv.writer(file)               # Name writer object
    writer.writerow(
        ["Accession", "Taxon ID", "Description", "Gene", "Sequence Length", "Date Late Modified", "Date Collected", "Domain", "Kingdom", "Superphylum", "Phylum",
         "Subphylum", "Class", "Subclass", "Infraclass", "Superorder", "Order", "Suborder", "Superfamily", "Family",
         "Subfamily", "Tribe", "Species", "Country", "Region", "Lat/Long", "Ref1 Author", "Ref1 Title", "Ref1 Journal", "Ref2 Author",
         "Ref2 Title", "Ref2 Journal", "Ref3 Author", "Ref3 Title", "Ref3 Journal"])  # Write column names


# Write row of metadata file
def writecsv(x):                                            # x = genbank record output
    row = [rec["gbid"], rec["txid"], rec["description"], rec["gene"], rec["length"], rec["rec date"], rec["c date"]]
    #while len(rec["taxonomy"]) <= 14:
        #rec[6].append("")
    rec["taxonomy"].extend([""] * (15 - len(rec["taxonomy"])))
    row.extend(rec["taxonomy"])
    row.append(rec["spec"])
    row.append(rec["country"])
    row.append(rec["region"])
    row.append(rec["latlon"])
    row.extend(rec["refs"])

    with open("metadata.csv", "a") as file:
        writer = csv.writer(file)
        writer.writerow(row)


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



genes = {"12S": ["12S", "12S RIBOSOMAL RNA", "12S RRNA"],
         "16S": ["16S", "16S RIBOSOMAL RNA", "16S RRNA"],
         "18S": ["18S", "18S RIBOSOMAL RNA", "18S RRNA", "18S SMALL SUBUNIT RIBOSOMAL RNA"],
         "EF1A": ["EF1-ALPHA", "EF1A", "ELONGATION FACTOR 1 ALPHA", "ELONGATION FACTOR 1-ALPHA"],
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
         "ND6": ['NAD6', 'NSD6', 'NADH6', 'NADH DEHYDROGENASE SUBUNIT VI', 'NADH DEHYDROGENASE SUBUNIT 6', 'NADH DESHYDROGENASE SUBUNIT 6', 'NAD6-0', 'ND6']}


# To use cli gene option, need to search entrez with list of all name variants.

#if args.gene:
#    geneslist = args.gene.split(",")
#    inputgenes = " OR ".join(geneslist)


unrecgenes = set()

# Set accepted genes and minimum sequence lengths
#min = {"ATP6": 500, "ATP8": 100, "COX1": 500, "COX2": 500, "COX3": 500, "CYTB": 500, "ND1": 500, "ND2": 500, "ND3": 300, "ND4": 500, "ND4L": 200, "ND5": 500, "ND6": 400}

x = 0  # Count taxids

accs = []
file = open("test.txt")
lines = file.readlines()
for line in lines:
    line.strip()
    accs.append(line)

accstr = ",".join(accs)                                           # Join into string for efetch
handle = Entrez.efetch(db="nucleotide", id=accstr, rettype="gb", retmode="text")  # Get GenBanks
record = SeqIO.parse(handle, "gb")
sequences = []
for rec in record:
    if "Dytiscidae" not in rec.annotations["taxonomy"]:
        continue
    for feature in rec.features:
        type = feature.type
        if type not in ('CDS', 'rRNA', 'mRNA'):
            continue  # skip the rest of the current iteration of this loop
        name = get_feat_name(feature)                       # Find gene name
        stdname = ""
        for k, v in genes.items():
            if name in v:
                stdname = k
        if stdname == "":
            unrecgenes.add(name)
            continue
        else:
            seq = feature.extract(rec.seq)
            if seq in sequences:
                continue
            else:
                sequences.append(seq)
                output = get_output(rec)
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
            break




print(f"\n{str(x)} gene records saved to species dict")

print("\nUnrecognised Genes")
print(unrecgenes)


# Set record length as 0, iterate through records and replace whenever another sequence is longer.
def findmax(x):
    max = x[0]["length"]
    maxrec = x[0]
    for record in x:
        if record["length"] > max:
            max = record["length"]
            maxrec = record
    return maxrec



# Dict for longest sequences, key is gene stdname, value is list of records
longest = {}
for tax, stdname in species.items():
    print(tax)
    for gene, records in stdname.items():
        print(gene)
        chosen = findmax(records)
        if gene in longest:
            longest[gene].append(chosen)
            print("add1")
        else:
            longest[gene] = [chosen]
            print("add2")

# Save each gene list to separate fasta file
# output = 0 gene, 1 GBID, 2 TXID, 3 description, 4 species, 5 date, 6 taxonomy(15 levels), 7 feature type,
#           8 sequence length, 9 sequence, 10 country, 11 region, 12 references
for gene, records in longest.items():
    file = open(f"{gene}.fasta", "w")
    for rec in records:
        #print(rec)
        writecsv(rec)
        file.write(f">{rec['gbid']}\n{rec['seq']}\n")
print("CSV and fastas written to file.")



