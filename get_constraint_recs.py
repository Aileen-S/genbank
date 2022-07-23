

#python3 get_other_genes.py -e mixedupvoyage@gmail.com -t Eretes


import argparse
import csv
from Bio import Entrez
from Bio import SeqIO
import textwrap as _textwrap
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation

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
parser = argparse.ArgumentParser(description="Search GenBank, retrieve gene sequences and save as fasta. Set up for 12S, 16S, 18S, 28S, AK, AS, CAD, EF1A, H3 and Wg", formatter_class=MultilineFormatter)
parser.add_argument("-t", "--taxon", type=str, help="Taxon of interest")
#parser.add_argument("-g", "--gene", type=str, help="Gene(s) of interest. Format: gene1,gene2,gene3")
parser.add_argument("-e", "--email", type=str, help="Your email registered with NCBI")


# Start the actual script

args = parser.parse_args()         # Process input args from command line
#args = argparse.Namespace(taxon='Amphizoidae', mpc=True, email='aileen.scott@nhm.ac.uk', nuclear=False) # This is how I step through the script interactively
#Namespace(taxon='Eretes', mpc=True, email='aileen.scott@nhm.ac.uk', nuclear=False)

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
species = {}
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
        z = 0
        gbid = rec.name
        if args.taxon not in rec.annotations["taxonomy"]:
            continue
        # Get record output data
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
        output = {"gbid": rec.name,
                  "txid": tax,
                  "description": rec.description,
                  "spec": rec.annotations["organism"],
                  "rec date": rec.annotations["date"],
                  "c date": c_date,
                  "taxonomy": rec.annotations["taxonomy"][0:15],
                  "country": country,
                  "region": region,
                  "latlon": latlon,
                  "refs": refs}
        #Get sequence output data
        for feature in rec.features:
            #sequence = rec[feature.location.start:feature.location.end]
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
                z += 1
            if "genes" in output:
                output["genes"][stdname] = {"gene": stdname,
                                            "type": type,
                                            "length": len(seq),
                                            "seq": seq}
            else:
                output["genes"] = {stdname: {"gene": stdname,
                                             "type": type,
                                             "length": len(seq),
                                             "seq": seq}}          #Temporarily changed
        if z >= 3:
            # only save recs with at least 3 genes
            if tax in species:
                if gbid in species[tax]:
                    species[tax][gbid].append(output)
                    x += 1
                else:
                    species[tax][gbid] = output
                    x += 1
            else:
                species[tax] = {gbid: output}
                x += 1
            break
print(f"\n{str(x)} gene records saved to species dict")

print("\nUnrecognised Genes")
print(unrecgenes)

# Write CSV metadata file
with open("metadata.csv", "w") as file:     # Open output file
    writer = csv.writer(file)               # Name writer object
    writer.writerow(
        ["Accession", "Taxon ID", "Description",
         '12S Length', '16S Length', '18S Length', 'EF1A Length', 'H3 Length', 'Wg Length', 'ATP6 Length',
         'ATP8 Length', 'COX1 Length', 'COX2 Length', 'COX3 Length', 'CYTB Length', 'ND1 Length', 'ND2 Length',
         'ND3 Length', 'ND4 Length', 'ND4L Length', 'ND5 Length', 'ND6 Length'
         "Date Late Modified", "Date Collected", "Domain", "Kingdom", "Superphylum", "Phylum",
         "Subphylum", "Class", "Subclass", "Infraclass", "Superorder", "Order", "Suborder", "Superfamily", "Family",
         "Subfamily", "Tribe", "Species", "Country", "Region", "Lat/Long", "Ref1 Author", "Ref1 Title", "Ref1 Journal",
         "Ref2 Author", "Ref2 Title", "Ref2 Journal", "Ref3 Author", "Ref3 Title", "Ref3 Journal"])

genes = ['12S', '16S', '18S', 'EF1A', 'H3', 'Wg', 'ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']
sequencedict = {}
file = open("metadata.csv", "a")
writer = csv.writer(file)
for rec in species.values():
    for gbid, output in rec.items():
        row = [output["gbid"], output["txid"], output["description"]]
        for gene in genes:
            if gene in output["genes"]:
                row.append(output["genes"][gene]["length"])
            else:
                row.append("")
    # Also include gene count column
        row.append(output["rec date"])
        row.append(output["c date"])
        output["taxonomy"].extend([""] * (15 - len(output["taxonomy"])))
        row.extend(output["taxonomy"])
        row.append(output["spec"])
        row.append(output["country"])
        row.append(output["region"])
        row.append(output["latlon"])
        row.extend(output["refs"])
        writer.writerow(row)
        for generec in output["genes"].values():            #write sequence dict for fastas
            if generec["gene"] in sequencedict:
                sequencedict[generec["gene"]][gbid] = generec["seq"]
            else:
                sequencedict[generec["gene"]] = {gbid: generec["seq"]}

for gene, gbids in sequencedict.items():
    file = open(f"{gene}.fasta", "w")
    for gbid, seq in gbids.items():
        file.write(f">{gbid} {gene}\n{seq}\n")
print("CSV and fastas written to file.")



