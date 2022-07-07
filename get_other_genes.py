

#python3 get_other_genes.py -e mixedupvoyage@gmail.com -t Eretes


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


# Argument parser
parser = argparse.ArgumentParser(description="Search GenBank, retrieve gene sequences and save as fasta. Set up for 12S, 16S, 18S, 28S, AK, AS, CAD, EF1A, H3 and Wg", formatter_class=MultilineFormatter)
parser.add_argument("-t", "--taxon", type=str, help="Taxon of interest")
#parser.add_argument("-g", "--gene", type=str, help="Gene(s) of interest. Format: gene1,gene2,gene3")
parser.add_argument("-e", "--email", type=str, help="Your email registered with NCBI")


# Start the actual script

args = parser.parse_args()         # Process input args from command line
#args = argparse.Namespace(taxon='Eretes', mpc=True, email='aileen.scott@nhm.ac.uk', nuclear=False) # This is how I step through the script interactively
#Namespace(taxon='Eretes', mpc=True, email='aileen.scott@nhm.ac.uk', nuclear=False)

genes = {"12S": ["12S RIBOSOMAL RNA", "12S RRNA"],
         "16S": ["16S RIBOSOMAL RNA", "16S RRNA"],
         "18S": ["18S RIBOSOMAL RNA", "18S RRNA", "18S SMALL SUBUNIT RIBOSOMAL RNA"],
         "28S": ["28S RIBOSOMAL RNA", "28S RRNA", "28S LARGE SUBUNIT RIBOSOMAL RNA"],
         "AK": ["AK", "ARGININE KINASE", "ARGK", "ARGKIN", "ARGS", "ARK"],
         "AS": [],
         "CAD": ["CAD", "CAD FRAGMENT 1"],
         "EF1A": ["EF1-ALPHA", "EF1A", "ELONGATION FACTOR 1 ALPHA", "ELONGATION FACTOR 1-ALPH"],
         "H3": ["H3"],
         "Wg": ["WG", "WINGLESS", "WNG", "WNT"]}


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
        for feature in rec.features:
            type = feature.type
            if type not in ('CDS', 'rRNA', 'mRNA'):
                continue  # skip the rest of the current iteration of this loop
            name = get_feat_name(feature)                       # Find gene name
            for k, v in genes.items():
                if name in v:
                    stdname = k
                    #if args.gene:
                        #if name not in geneslist:
                            #.add(name)
                            #continue
                    sequence = rec[feature.location.start:feature.location.end]
                    #if stdname in min:
                        #if len(sequence) < min[stdname]:
                            #continue
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
                    #print(type(country))
                    seq = bytes(sequence.seq)   # Tried changing to bytes
                    if seq == "None":           # and skip features with no sequence.
                        continue                # Still getting "Bio.Seq.UndefinedSequenceError"
                    output = {"gene" : stdname,
                              "gbid" : rec.name,
                              "txid" : tax,
                              "description" : rec.description,
                              "spec" : rec.annotations["organism"],
                              "rec date" : rec.annotations["date"],
                              "c date" : c_date,
                              "taxonomy" : rec.annotations["taxonomy"][0:15],
                              "type" : type,
                              "length" : len(sequence),
                              "seq" : bytes.decode(seq),
                              "country" : country,
                              "region" : region,
                              "latlon" : latlon,
                              "refs" : refs}
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

print(f"\n{str(x)} gene records saved to species dict")
#print(species)

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
    for gene, records in stdname.items():
        chosen = findmax(records)
        if gene in longest:
            longest[gene].append(chosen)
        else:
            longest[gene] = [chosen]


# Save each gene list to separate fasta file
# output = 0 gene, 1 GBID, 2 TXID, 3 description, 4 species, 5 date, 6 taxonomy(15 levels), 7 feature type,
#           8 sequence length, 9 sequence, 10 country, 11 region, 12 references
for gene, records in longest.items():
    file = open(f"{gene}.fa", "w")
    for rec in records:
        #print(rec)
        writecsv(rec)
        file.write(f">{rec['gbid']}\n{rec['seq']}\n")

