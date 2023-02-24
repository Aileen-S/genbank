
#python3 get_concat_recs.py -e mixedupvoyage@gmail.com -t Megadytes

#python3 get_concat_recs.py -e mixedupvoyage@gmail.com -f test.txt -r gbid -i both

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
# Add option to find only mito genes, or only selected genes.
parser = argparse.ArgumentParser(description="Search GenBank, retrieve COI sequences and save as fasta.")
parser.add_argument('-f', '--file', type=str, help="Input file with accession or taxon ID list")
parser.add_argument('-i', '--fasta_id', choices=['gbid', 'txid'], help="Choose identifier for output fastas. Default is both txid and gbid.")

args = parser.parse_args()         # Process input args from command line
#args = argparse.Namespace(taxon='Amphizoidae', mpc=True, email='aileen.scott@nhm.ac.uk', nuclear=False) # This is how I step through the script interactively


coi = ['CYTOCHROME C OXIDASE SUBUNIT 1', 'CYTOCHROME OXIDASE SUBUNIT I', 'CYTOCHROME C OXIDASE SUBUNIT I', 'COXI', 'CO1', 'COI', 'CYTOCHROME COXIDASE SUBUNIT I', 'CYTOCHROME OXIDASE SUBUNIT 1', 'CYTOCHROME OXYDASE SUBUNIT 1', 'COX1']

#CYTOCHROME C OXIDASE SUBUNIT 1 OR CYTOCHROME OXIDASE SUBUNIT I OR CYTOCHROME C OXIDASE SUBUNIT I OR COXI OR CO1 OR COI OR CYTOCHROME COXIDASE SUBUNIT I OR CYTOCHROME OXIDASE SUBUNIT 1 OR CYTOCHROME OXYDASE SUBUNIT 1 OR COX1

alltxids = open("alltxids.csv", 'w')
allwriter = csv.writer(alltxids)  # Name writer object
allwriter.writerow(["Taxon ID", "Suborder", "Superfamily", "Family", "Subfamily", "Tribe", 'Genus', 'Species', "Description"])

coitxids = open('coitxids.csv', 'w')
coiwriter = csv.writer(coitxids)  # Name writer object
coiwriter.writerow(["Taxon ID", "Suborder", "Superfamily", "Family", "Subfamily", "Tribe", 'Genus', 'Species', "Description"])

unrec_genes = set()

# Search through GBIDs
print(f'Searching {args.file} for COI sequences')


species = {}
x = 0  # Count taxids
r = 0  # Count records
with open(args.file) as file:      # Specify location and name file
    for rec in SeqIO.parse(file, "gb"):
        r += 1
        if r % 1000 == 0:
            print(f'{r} records searched')
        db_xref = rec.features[0].qualifiers["db_xref"]
        for ref in db_xref:
            if "taxon" in ref:  # Get NCBI taxon, rather than BOLD cross ref
                txid = "".join(filter(str.isdigit, ref))  # Extract numbers from NCBI taxon value
        spec = rec.annotations["organism"]
        specfasta = spec.replace(" ", "_")
        taxonomy = rec.annotations["taxonomy"][10:15]
        taxonomy.extend([""] * (5 - len(taxonomy)))
        fastatax = f"{taxonomy[2]}_{taxonomy[3]}_{taxonomy[4]}_{specfasta}"
        if "country" in rec.features[0].qualifiers:
            location = rec.features[0].qualifiers["country"][0]
            if ":" in location:
                country, region = location.split(":", 1)
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
        if "references" in rec.annotations:
            for ref in rec.annotations["references"]:
                refs.append(ref.authors)
                refs.append(ref.title)
                refs.append(ref.journal)
        row = [txid]
        row.extend(taxonomy)
        gen_spec = spec.split(' ')
        row.append(gen_spec[0])
        row.append(gen_spec[1])
        row.append(rec.description)
        allwriter.writerow(row)
        for feature in rec.features:
            type = feature.type
            if type not in ('CDS', 'rRNA'):
                continue  # skip the rest of the current iteration of this loop
            name = get_feat_name(feature)                       # Find gene name
            stdname = ''
            if name not in coi:
                unrec_genes.add(name)
                continue
            seq = feature.extract(rec.seq)
            if 657 <= len(seq) <= 658:
                continue
            if 'codon_start' in feature.qualifiers:
                frame = feature.qualifiers["codon_start"]
            else:
                frame = ''
                print(f"Reading frame missing from record {rec.name}, {stdname}.")
            output = {"gene": stdname,
                      "gbid": rec.name,
                      "txid": txid,
                      "description": rec.description,
                      "spec": rec.annotations["organism"],
                      "rec date": rec.annotations["date"],
                      "c date": c_date,
                      "taxonomy": taxonomy,
                      "fastatax": fastatax,
                      "type": type,
                      "length": len(seq),
                      "seq": seq,
                      "frame": frame,
                      "country": country,
                      "region": region,
                      "latlon": latlon,
                      "refs": refs}
            crow = [output["txid"]]
            crow.extend(output["taxonomy"])
            gen_spec = output['spec'].split(' ')
            crow.append(gen_spec[0])
            crow.append(gen_spec[1])
            crow.append(output['description'])
            coiwriter.writerow(crow)

            if txid in species:                              # If taxon ID in dict
                species[txid].append(output)    # Add gene info list to dict
                x += 1
            else:
                species[txid] = [output]      # Otherwise add to dict with new key
                x += 1
            #break

print(f'\nTotal {r} records searched')
print(f"{x} COI sequences found")
print(f'{len(species.keys())} taxon IDs found')

print("\nUnrecognised Genes")
print(f'{unrec_genes}\n')

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
longest = []
for tax, records in species.items():
    chosen = findmax(records)
    longest.append(chosen)
# Save each gene list to separate fasta file

# Write CSV metadata file
with open("metadata.csv", "w") as file:     # Open output file
    writer = csv.writer(file)               # Name writer object
    writer.writerow(
        ["Accession", "Taxon ID", "Species", 'COI Length',
         "Suborder", "Superfamily", "Family", "Subfamily", "Tribe", 'Genus', "Description", "Date Late Modified",
         "Date Collected", "Country", "Region", "Lat/Long", "Ref1 Author", "Ref1 Title", "Ref1 Journal", "Ref2 Author",
         "Ref2 Title", "Ref2 Journal", "Ref3 Author", "Ref3 Title", "Ref3 Journal"])

file = open(f"{args.input.split('.')[0]}_metadata.csv", "a")
writer = csv.writer(file)

for output in longest:
    row = [output["gbid"], output["txid"], output["spec"], output['length']]
    row.extend(output["taxonomy"])
    gen_spec = output['spec'].split(' ')
    genus = gen_spec[0]
    row.append(genus)
    row.append(output["description"])
    row.append(output["rec date"])
    row.append(output["c date"])
    row.append(output["country"])
    row.append(output["region"])
    row.append(output["latlon"])
    row.extend(output["refs"])
    writer.writerow(row)

file = open(f"{args.input.split('.')[0]}.fasta", "w")
x = 0
for rec in longest:
    if args.fasta_id:
        if args.fasta_id == 'txid':
            f_id = rec['txid']
        if args.fasta_id == 'gbid':
            f_id = rec['gbid']
    else:
        f_id = f"{rec['txid']}_{rec['gbid']}"
    file.write(f">{f_id}_{rec['fastatax']};frame={rec['frame'][0]}\n{rec['seq']}\n")
    x += 1

print(f'{x} records written to COI.fasta')
print("Metadata saved to metadata.csv")



