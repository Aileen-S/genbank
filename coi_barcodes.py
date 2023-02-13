
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
parser = argparse.ArgumentParser(description="Search GenBank, retrieve gene sequences and save as fasta.")
parser.add_argument("-t", "--taxon", type=str, help="Taxon of interest")
parser.add_argument('-f', '--file', type=str, help="Input file with accession or taxon ID list")
parser.add_argument('-r', '--ref', choices=['txid', 'gbid'], help="If using --file option, specify accessions or taxon IDs.")
parser.add_argument('-i', '--fasta_id', choices=['gbid', 'txid', 'both'], help="Choose identifiers for output fastas. Default is gbid.")
parser.add_argument("-e", "--email", type=str, help="Your email registered with NCBI")


args = parser.parse_args()         # Process input args from command line
#args = argparse.Namespace(taxon='Amphizoidae', mpc=True, email='aileen.scott@nhm.ac.uk', nuclear=False) # This is how I step through the script interactively


coi = ['CYTOCHROME C OXIDASE SUBUNIT 1', 'CYTOCHROME OXIDASE SUBUNIT I', 'CYTOCHROME C OXIDASE SUBUNIT I', 'COXI', 'CO1', 'COI', 'CYTOCHROME COXIDASE SUBUNIT I', 'CYTOCHROME OXIDASE SUBUNIT 1', 'CYTOCHROME OXYDASE SUBUNIT 1', 'COX1']


unrec_genes = set()
unrec_species = []
Entrez.email = args.email

if args.ref == 'gbid':
    accs = []
    file = open(args.file)
    lines = file.readlines()
    for line in lines:
        acc = line.strip()
        accs.append(acc)

else:
    if args.ref  == 'txid':
        taxids = []
        file = open(args.file)
        lines = file.readlines()
        for line in lines:
            taxid = line.strip()
            taxids.append(taxid)

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

    txfile = open('txids.txt', 'w')
    for txid in taxids:
        txfile.write(f'{txid}\n')
    print('Taxon IDs saved to txids.txt')

    y = 0  # Count records saved
    accs = []
    for tax in taxids:
        y += 1
        if y % 100 == 0:
            print(f"Downloading GenBank records for taxon IDs {y+1} to {y+100}" if (y+100) < len(taxids) else
                  f"Downloading GenBank records for taxon IDs {y+1} to {len(taxids)}")
        handle = Entrez.esearch(db="nucleotide", term=f"txid{tax}")       # Search for all records for each taxon id
        record = Entrez.read(handle)
        accs   = accs + record["IdList"]   # Get GBIDs

# Search through GBIDs
species = {}
x = 0  # Count taxids
accstr = ",".join(accs)                                           # Join into string for efetch
handle = Entrez.efetch(db="nucleotide", id=accstr, rettype="gb", retmode="text")  # Get GenBanks
record = SeqIO.parse(handle, "gb")
for rec in record:
    if args.taxon:
        if args.taxon not in rec.annotations["taxonomy"]:
            unrec_species.append(rec.name)
            continue
    db_xref = rec.features[0].qualifiers["db_xref"]
    for ref in db_xref:
        if "taxon" in ref:  # Get NCBI taxon, rather than BOLD cross ref
            txid = "".join(filter(str.isdigit, ref))  # Extract numbers from NCBI taxon value
    spec = rec.annotations["organism"]
    specfasta = spec.replace(" ", "_")
    taxonomy = rec.annotations["taxonomy"][10:15]
    taxonomy.extend([""] * (5 - len(taxonomy)))
    if taxonomy[4] == "Cybistrini":
        taxonomy[3] = "Cybistrinae"
    fastatax = f"{taxonomy[2]}_{taxonomy[3]}_{taxonomy[4]}_{specfasta}"

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
    if "references" in rec.annotations:
        for ref in rec.annotations["references"]:
            refs.append(ref.authors)
            refs.append(ref.title)
            refs.append(ref.journal)
    else:
        continue
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
        if 680 <= len(seq) <= 1000:
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
        if txid in species:                              # If taxon ID in dict
                species[txid].append(output)    # Add gene info list to dict
                x += 1
        else:
            species[txid] = [output]      # Otherwise add to dict with new key
            x += 1
        #break

print(f"\n{x} gene records saved to species dict")

print("\nUnrecognised Genes")
print(f'{unrec_genes}\n')
#print("\nUnrecognised Species")
#print(unrec_species)


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

file = open("metadata.csv", "a")
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

file = open(f"COI.fasta", "w")
for rec in longest:
    x = 0
    if args.fasta_id:
        if args.fasta_id == 'txid':
            f_id = rec['txid']
        if args.fasta_id == 'both':
            f_id = f"{rec['txid']}_{rec['gbid']}"
    else:
        f_id = rec['gbid']
    file.write(f">{f_id}_{rec['fastatax']};frame={rec['frame'][0]}\n{rec['seq']}\n")
    x += 1

print(f'{x} records written to COI.fasta')
print("Metadata saved to metadata.csv")



