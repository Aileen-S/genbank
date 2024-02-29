
#python3 get_concat_recs.py -e mixedupvoyage@gmail.com -t Megadytes

import argparse
import csv
from Bio import SeqIO


# Function definitions

def get_feat_name(feat):
    featname = "unknown"
    nametags = ['gene', 'product', 'label', 'standard_name', 'note']  # Search these keys for gene name
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


# Argument parser
# Add option to find only mito genes, or only selected genes.
parser = argparse.ArgumentParser(description="Search GenBank file, retrieve gene sequences and save as fasta.")
parser.add_argument("-t", "--taxon", type=str, help="Taxon of interest")
parser.add_argument('-g', '--gb_file', type=str, help="Input genbank format file")
parser.add_argument('-a', '--accs', type=str, help="Input file with list of genbank accession numbers")
parser.add_argument('-x', '--txid', type=str, help="Input file with list of genbank taxon IDs")

parser.add_argument('-m', '--mito', action='store_true', help='Save only mitochondrial protein-coding genes')
parser.add_argument('-c', '--coi', action='store_true', help='Save only COX1')

parser.add_argument('-l', '--longest', action='store_true', help='Save only longest sequences per gene per taxon ID')

parser.add_argument('-i', '--fasta_id', choices=['gbid', 'txid', 'both'], help="Choose identifiers for output fastas. Default is gbid.")
parser.add_argument('-s', '--skip', type=str, help="File with list of GBIDs to avoid")


args = parser.parse_args()         # Process input args from command line
#args = argparse.Namespace(taxon='Amphizoidae', mpc=True, email='aileen.scott@nhm.ac.uk', nuclear=False) # This is how I step through the script interactively

genes = {"12S": ["12S", "12S RIBOSOMAL RNA", "12S RRNA", "RRNS", "SSU"],
         "16S": ["16S", "16S RIBOSOMAL RNA", "16S RRNA", "RRNL", "LARGESUBUNITRIBOSOMALRNA", "LSU", "SSRNL"],
         "ATP6": ['ATP SYNTHASE F0 SUBUNIT 6', 'APT6', 'ATP SYNTHASE A0 SUBUNIT 6', 'ATP SYNTHASE SUBUNIT 6', 'ATP SYNTHASE FO SUBUNIT 6', 'ATPASE6', 'ATPASE SUBUNIT 6', 'ATP6'],
         "ATP8": ['ATP SYNTHASE F0 SUBUNIT 8', 'APT8', 'ATP SYNTHASE A0 SUBUNIT 8', 'ATP SYNTHASE SUBUNIT 8', 'ATP SYNTHASE FO SUBUNIT 8', 'ATPASE8', 'ATPASE SUBUNIT 8', 'ATP8'],
         "COX1": ['CYTOCHROME C OXIDASE SUBUNIT 1', 'CYTOCHROME OXIDASE SUBUNIT I', 'CYTOCHROME C OXIDASE SUBUNIT I', 'COXI', 'CO1', 'COI', 'CYTOCHROME COXIDASE SUBUNIT I', 'CYTOCHROME OXIDASE SUBUNIT 1', 'CYTOCHROME OXYDASE SUBUNIT 1', 'COX1', 'CYTOCHROME OXIDASE I', 'CYTOCHROME OXIDASE C SUBUNIT I', 'COX 1'],
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

misc = ["similar to cytochrome oxidase subunit I", "similar to cytochrome oxidase subunit 1", "similar to cytochrome c oxidase subunit I; COI", "similar to cytochrome c oxidase subunit I", "sequence contains partial cox1, tRNA-Leu and partial cox2 genes", "cox1", "coding region not determined; cytochrome oxidase subunit I", "similar to cytochrome oxidase subunit I; COI", "similar to cytochrome c oxidase", "similar to cytochrome oxidase subunit 1; 3' barcoding region; LCO-HCO"]

suborders = ['Adephaga', 'Polyphaga', 'Myxophaga', 'Archostemata']

unrec_genes = set()
unrec_species = []
misc_ids = []
misc_sequences = []

if args.accs:
    accs = []
    file = open(args.accs)
    lines = file.readlines()
    for line in lines:
        acc = line.strip()
        accs.append(acc)
    print(f'{len(accs)} IDs found in {args.accs}')

if args.skip:
    skip = []
    file = open(args.skip)
    lines = file.readlines()
    for line in lines:
        rec = line.strip()
        skip.append(rec)
    print(f'{len(skip)} IDs to avoid found in {args.skip}')


if args.txid:
    txids = []
    file = open(args.txid)
    lines = file.readlines()
    for line in lines:
        txid = line.strip()
        txids.append(txid)
    print(txids)

# Search through GBIDs
count = open('count.csv', 'w')
species = {}
sequences = []
nohits = []
other_type = set()
misc_feature = set()
x = 0  # Count taxids
with open(args.gb_file) as file:
    record = SeqIO.parse(file, "gb")
    for rec in record:
        if args.accs:
            if rec.name not in accs:
                continue
        if args.taxon:
            if args.taxon not in rec.annotations["taxonomy"]:
                unrec_species.append(rec.name)
                continue
        try:
            db_xref = rec.features[0].qualifiers["db_xref"]
        except KeyError:
            continue
        bold = ''
        for ref in db_xref:
            if "taxon" in ref:  # Get NCBI taxon, rather than BOLD cross ref
                txid = "".join(filter(str.isdigit, ref))  # Extract numbers from NCBI taxon value
            if "BOLD" in ref:
                bold = (ref.split(":")[1]).split('.')[0]
        if args.txid:
            if txid not in txids:
                continue
        if args.skip:
            if txid in skip:
                continue
        spec = rec.annotations["organism"]
        # Replace the following characters: > < . ( ) ; : ' ,
        spec = spec.replace(">", "_").replace("<", "_").replace(".", "").replace('(', '_')\
            .replace(')', '_').replace(';', '_').replace(':', '_').replace("'", "").replace(',', '')
        specfasta = spec.replace(" ", "_")
        taxonomy = ['', '', '', '', '']
        for tax in rec.annotations["taxonomy"]:
            if tax in suborders: taxonomy[0] = tax
            if tax.endswith('oidea'): taxonomy[1] = tax
            if tax.endswith('idae'): taxonomy[2] = tax
            if tax.endswith('inae'): taxonomy[3] = tax
            if tax.endswith('ini'): taxonomy[4] = tax
        taxonomy.append(spec.split(' ')[0])
        fastatax = f"{taxonomy[2]}_{taxonomy[3]}_{taxonomy[4]}_{specfasta}"
        if args.fasta_id:
            if args.fasta_id == 'txid':
                f_id = txid
            if args.fasta_id == 'both':
                f_id = f"{txid}/_{rec.name}/"
            else:
                f_id = f"{rec.name}/"
        else:
            f_id = f"{rec.name}/"
        fasta_id = f'{f_id}_{fastatax}'
        if "country" in rec.features[0].qualifiers:
            location = rec.features[0].qualifiers["country"][0]
            if ":" in location:
                country, region = location.split(":",1)
            else:
                country = location
                region = ""
        else:
            country = ""
            region = ""
        if "lat_lon" in rec.features[0].qualifiers:
            latlon = rec.features[0].qualifiers["lat_lon"][0]
            ll_list = latlon.split(" ")
            if ll_list[1] == "N":
                lat = ll_list[0]
            else:
                lat = "-" + ll_list[0]
            if ll_list[3] == "E":
                long = ll_list[2]
            else:
                long = "-" + ll_list[2]
        else:
            latlon = ""
            lat = ""
            long = ""
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
            references = ''
        g = 0
        for feature in rec.features:
            type = feature.type
            if type not in ('CDS', 'rRNA'):
                continue
                #if type == 'misc_feature':
                #    try:
                #        if feature.qualifiers['note'][0] in misc:
                #            seq = feature.extract(rec.seq)
                #            misc_ids.append(rec.name)
                #            stdname = 'COX1'
                #            rec.name = f'misc_{rec.name}'
                #            frame = ['']
                #    except KeyError:
                #        continue
                #else:
                #    other_type.add(type)
                #    continue
            else:
                name = get_feat_name(feature)                       # Find gene name
                stdname = ""
                for k, v in genes.items():
                    if name in v:
                        stdname = k
                        g += 1
                if stdname == '':
                    unrec_genes.add(name)
                    continue
                if args.mito:
                    if stdname not in mito:
                        continue
                if args.coi:
                    if stdname != 'COX1':
                        continue
                if stdname in cds:
                    if 'codon_start' in feature.qualifiers:
                        frame = feature.qualifiers["codon_start"]
                    else:
                        frame = ['']
                        print(f"Reading frame missing from record {rec.name}, {stdname}.")
                else:
                    frame = ''
            seq = ''
            seq = feature.extract(rec.seq)
            sequences.append(seq)
            output = {"gene": stdname,
                      "gbid": rec.name,
                      "txid": txid,
                      "bold": bold,
                      "description": rec.description,
                      "spec": spec,
                      "rec date": rec.annotations["date"],
                      "c date": c_date,
                      "taxonomy": taxonomy,
                      "fasta_id": fasta_id,
                      "type": type,
                      "length": len(seq),
                      "seq": seq,
                      "frame": frame,
                      "country": country,
                      "region": region,
                      "latlon": latlon,
                      "lat": lat,
                      "long": long,
                      "refs": refs}
            #if rec.name in misc_ids:
            #    misc_sequences.append(output)
            #else:
            if txid in species:                              # If taxon ID in dict
                if stdname in species[txid]:                 # If gene in dict for that taxon ID
                    species[txid][stdname].append(output)    # Add gene info list to dict
                    x += 1
                else:
                    species[txid][stdname] = [output]      # Otherwise add to dict with new key
                    x += 1
            else:
                species[txid] = {stdname: [output]}      # Otherwise add to dict with new key
                x += 1
        if g == 0:
            nohits.append(rec.name)
        else:
            count.write(f'{rec.name},{g}\n')



print(f"\n{x} sequences found for requested genes\n")


# Set record length as 0, iterate through records and replace whenever another sequence is longer.
def findmax(x):
    max = x[0]["length"]
    maxrec = x[0]
    for record in x:
        if record["length"] > max:
            max = record["length"]
            maxrec = record
    return maxrec


# Write CSV metadata file
with open("metadata.csv", "w") as file:     # Open output file
    writer = csv.writer(file)               # Name writer object
    writer.writerow(
        ["FastaID", "Accession", "TXID", 'BOLD', "Species", 'Gene', 'Length',
         "Suborder", "Superfamily", "Family", "Subfamily", "Tribe", 'Genus', "Description", "Date Late Modified",
         "Date Collected", "Country", "Region", "Lat/Long", "Lat", "Long", "Ref1 Author", "Ref1 Title", "Ref1 Journal",
         "Ref2 Author", "Ref2 Title", "Ref2 Journal", "Ref3 Author", "Ref3 Title", "Ref3 Journal"])

gen = ['ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']


# Dict for longest sequences, key is gene stdname, value is list of records
longest = {}
for tax, stdname in species.items():
    for gene, records in stdname.items():
        # Get longest sequence for each taxon ID, for each gene
        if args.longest:
            chosen = findmax(records)
            if gene in longest:
                longest[gene].append(chosen)
            else:
                longest[gene] = [chosen]
        # Keep all sequences
        else:
            for record in records:
                if gene in longest:
                    longest[gene].append(record)
                else:
                    longest[gene] = [record]

# Save each gene list to separate fasta file

file = open("metadata.csv", "a")
writer = csv.writer(file)
for gene, records in longest.items():
    for output in records:
        row = [output["fasta_id"], output["gbid"], output["txid"], output["bold"], output["spec"], output['gene'], output['length']]

    # Also include gene count column
        row.extend(output["taxonomy"])
        row.append(output["description"])
        row.append(output["rec date"])
        row.append(output["c date"])
        row.append(output["country"])
        row.append(output["region"])
        row.append(output["latlon"])
        row.append(output["lat"])
        row.append(output["long"])
        row.extend(output["refs"])
        writer.writerow(row)

for output in misc_sequences:
    row = [output["gbid"], output["txid"], output["bold"], output["spec"], 'yes']
    row.extend(['', '', output['length'], '', '', '', '', '', '', '', '', '', ''])
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


for gene, records in longest.items():
    file = open(f"{gene}.fasta", "w")
    x = 0
    for rec in records:
        if gene in cds:
            file.write(f">{rec['fasta_id']};frame={rec['frame'][0]}\n{rec['seq']}\n")
            x += 1
        else:
            file.write(f">{rec['fasta_id']}\n{rec['seq']}\n")
            x += 1
    print(f'{x} records written to {gene}.fasta')

#file = open('misc.fasta', 'w')
#x = 0
#for misc in misc_sequences:
#    if args.fasta_id:
#        if args.fasta_id == 'txid':
#            f_id = f"{misc['txid']}/"
#        if args.fasta_id == 'both':
#            f_id = f"{misc['txid']}/_{misc['gbid']}/"
#        else:
#            f_id = f"{misc['gbid']}/"
#    else:
#        f_id = f"{misc['gbid']}/"
#    file.write(f">{f_id}_{misc['fastatax']}\n{misc['seq']}\n")
#    x += 1
#print(f'{x} records written to misc.fasta')


print("Metadata saved to metadata.csv")

if args.accs:
    if len(nohits) > 0:
        print(f'\nNo requested genes found in the following records: {nohits}')

print("\nUnrecognised Genes")
print(f'{unrec_genes}\n')
print('Misc Features')
print(misc_feature)
print("Other Feature Types")
print(other_type)