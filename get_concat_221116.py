
#python3 get_concat_recs.py -e mixedupvoyage@gmail.com -t Megadytes

# Make dict of subgenus names: new colum for genus, renamed from subgenus dict if necessary.

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
parser.add_argument('-i', '--fasta_id', choices=['gbid', 'txid', 'both'], help="Choose identifiers for output fastas. Default is gbid.")
parser.add_argument('-b', '--both', action="store_true", help="Print taxon ID and accession in output fastas.")
parser.add_argument("-e", "--email", type=str, help="Your email registered with NCBI")


args = parser.parse_args()         # Process input args from command line
#args = argparse.Namespace(taxon='Amphizoidae', mpc=True, email='aileen.scott@nhm.ac.uk', nuclear=False) # This is how I step through the script interactively

genes = {"12S": ["12S", "12S RIBOSOMAL RNA", "12S RRNA"],
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
         "28S": ["28S RIBOSOMAL RNA", "28S RRNA", "28S LARGE SUBUNIT RIBOSOMAL RNA"],
         "AK": ["AK", "ARGININE KINASE", "ARGK", "ARGKIN", "ARGS", "ARK"],
         "CAD": ["CAD", "CAD FRAGMENT 1", "CARBAMOYLPHOSPHATE SYNTHETASE"],
         "EF1A": ["EF1-ALPHA", "EF1A", "ELONGATION FACTOR 1 ALPHA", "ELONGATION FACTOR 1-ALPHA"],
         "H3": ["H3", "HISTONE 3", "HISTONE H3", "HIS3"],
         "RNApol": ["RNA POL II", "RNA POL2", "RNA POLYMERASE II LARGE SUBUNIT"],
         "Wg": ["WG", "WINGLESS", "WNG", "WNT", "WNT1", "WNT-4"]}

mito = ['ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']
nuc = ['AK', 'CAD', 'EF1A', 'H3', 'RNApol', 'Wg']
rna = ['12S', '16S', '18S', '28S']
cds = ['ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6', 'AK', 'CAD', 'EF1A', 'H3', 'RNApol', 'Wg']


subtribes = {'Deronectina': ['Amurodytes', 'Boreonectes', 'Clarkhydrus', 'Deronectes', 'Deuteronectes', 'Hornectes', 'Iberonectes',
                             'Larsonectes', 'Leconectes', 'Mystonectes', 'Nebrioporus', 'Nectoboreus', 'Nectomimus', 'Nectoporus',
                             'Neonectes', 'Oreodytes', 'Scarodytes', 'Stictotarsus', 'Trichonectes', 'Zaitzevhydrus'],
             'Hydroporina': ['Haideoporus', 'Heterosternuta', 'Hydrocolus', 'Hydroporus', 'Neoporus', 'Sanfilippodytes'],
             'Siettitiina': ['Ereboporus', 'Etruscodytes', 'Graptodytes', 'Iberoporus', 'Lioporeus', 'Metaporus', 'Porhydrus',
                             'Psychopomporus', 'Rhithrodytes', 'Siettitia', 'Stictonectes', 'Stygoporus'],
             'Sternopriscina': ['Antiporus', 'Barretthydrus', 'Brancuporus', 'Carabhydrus', 'Chostonectes', 'Megaporus',
                                'Necterosoma', 'Paroster', 'Sekaliporus', 'Sternopriscus', 'Tiporus']}

unrec_genes = set()
unrec_species = []
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
sequences = []
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
        if args.taxon not in rec.annotations["taxonomy"]:
            unrec_species.append(rec.name)
            continue
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
        for ref in rec.annotations["references"]:
            refs.append(ref.authors)
            refs.append(ref.title)
            refs.append(ref.journal)

        for feature in rec.features:
            type = feature.type
            if type not in ('CDS', 'rRNA'):
                continue  # skip the rest of the current iteration of this loop
            name = get_feat_name(feature)                       # Find gene name
            stdname = ""
            for k, v in genes.items():
                if name in v:
                    stdname = k
            if stdname == '':
                continue
            if stdname == "":
                unrec_genes.add(name)
                continue
            if stdname in cds:
                if 'codon_start' in feature.qualifiers:
                    frame = feature.qualifiers["codon_start"]
                else:
                    frame = ''
                    print(f"Reading frame missing from record {rec.name}, {stdname}.")
            else:
                frame = ''
            seq = feature.extract(rec.seq)
            sequences.append(seq)
            output = {"gene": stdname,
                      "gbid": rec.name,
                      "txid": tax,
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
            #break


print(f"\n{x} gene records saved to species dict")

print("\nUnrecognised Genes")
print(unrec_genes)
print("\nUnrecognised Species")
print(unrec_species)


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

# Write CSV metadata file
with open("metadata.csv", "w") as file:     # Open output file
    writer = csv.writer(file)               # Name writer object
    writer.writerow(
        ["Accession", "Taxon ID", "Species", '18S', "28S", "AK", "CAD", 'EF1A', 'H3', 'RNApol', 'Wg',
         '12S', '16S', 'ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6',
         "Suborder", "Superfamily", "Family", "Subfamily", "Tribe", 'Subtribe', 'Genus', "Description", "Date Late Modified",
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


file = open("metadata.csv", "a")
writer = csv.writer(file)
for gene, records in longest.items():
    for output in records:
        row = [output["gbid"], output["txid"], output["spec"]]
        for g in gen:
            if g == output['gene']:
                row.append(output['length'])
            else:
                row.append("")
    # Also include gene count column

        row.extend(output["taxonomy"])
        gen_spec = output['spec'].split(' ')
        genus = gen_spec[0]
        for k, v in subgenus.items():
            if genus in v:
                genus = k
        for k, v in subtribes.items():
            if genus in v:
                subtribe = k
            else:
                subtribe = ''
        row.append(subtribe)
        row.append(genus)
        row.append(output["description"])
        row.append(output["rec date"])
        row.append(output["c date"])
        row.append(output["country"])
        row.append(output["region"])
        row.append(output["latlon"])
        row.extend(output["refs"])
        writer.writerow(row)



for gene, records in longest.items():
    file = open(f"{gene}.fasta", "w")
    x = 0
    for rec in records:
        if args.fasta_id:
            if args.fasta_id == 'txid':
                f_id = rec['txid']
            if args.fasta_id == 'both':
                f_id = f"{rec['txid']}_{rec['gbid']}"
            else:
                f_id = rec['gbid']
        else:
            f_id = rec['gbid']
        if gene in cds:
            file.write(f">{f_id}_{rec['fastatax']};frame={rec['frame'][0]}\n{rec['seq']}\n")
            x += 1
        else:
            file.write(f">{f_id}_{rec['fastatax']}\n{rec['seq']}\n")
            x += 1
    print(f'{x} records written to {gene}.fasta')


print("Metadata saved to metadata.csv")


