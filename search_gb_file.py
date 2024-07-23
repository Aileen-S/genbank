
#python3 get_concat_recs.py -e mixedupvoyage@gmail.com -t Megadytes

import argparse, argcomplete
import csv
from Bio import SeqIO
from collections import Counter
import re

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

def genbank_metadata(rec):
    # NCBI taxon ID
    db_xref = rec.features[0].qualifiers.get("db_xref", [])
    txid = ""
    for ref in db_xref:
        if "taxon" in ref:  # Get NCBI taxon, rather than BOLD cross ref
            txid = "".join(filter(str.isdigit, ref))  # Extract numbers from NCBI taxon value

    # Taxonomy
    # Replace the following characters: > < . ( ) ; : ' ,
    spec = re.sub(r"[><.();:'\"]", "", rec.annotations["organism"]).replace(",", "")
    spec_parts = [part for part in spec.split(" ") if not re.search(r'\d', part) and not part.isupper()]
    spec = " ".join(spec_parts)
    specfasta = spec.replace(" ", "_")

    taxonomy = ['', '', '', '', '', '']
    for tax in rec.annotations["taxonomy"]:
        if tax in suborders: taxonomy[0] = tax
        if tax.endswith('formia'): taxonomy[1] = tax
        if tax.endswith('oidea'): taxonomy[2] = tax
        if tax.endswith('idae'): taxonomy[3] = tax
        if tax.endswith('inae'): taxonomy[4] = tax
        if tax.endswith('ini'): taxonomy[5] = tax
    #taxonomy.append(spec.split(' ')[0])
    #fastatax = f"{txid}_{taxonomy[2]}_{taxonomy[3]}_{taxonomy[4]}_{specfasta}"

    # Location
    if "country" in rec.features[0].qualifiers:
        location = rec.features[0].qualifiers["country"][0]
        if ":" in location:
            country, region = location.split(":", 1)
        else:
            country = location
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

    # References
    refs = []
    if "references" in rec.annotations:
        first = rec.annotations['references'][0]
        refs.append(first.authors)
        refs.append(first.title)
        refs.append(first.journal)
    output = {"gbid": rec.name,
              "txid": txid,
              "description": rec.description,
              "spec_id": rec.annotations["organism"],
              "spec": spec,
              "date": rec.annotations["date"],
              "taxonomy": taxonomy,
              "country": country,
              "lat": lat,
              "long": long,
              "refs": refs,
              "row": [txid, rec.name, '', '', ''] + taxonomy + [rec.annotations["organism"], country, lat, long] + refs}
    return output



# Argument parser
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

argcomplete.autocomplete(parser)
args = parser.parse_args()

genes = {"12S": ["12S", "12S RIBOSOMAL RNA", "12S RRNA", "RRNS", "SSU", "RRN12", "S-RRNA", "12S SMALL SUBUNIT RIBOSOMAL RNA", "SMALL SUBUNIT RIBOSOMAL RNA"],
         "16S": ["16S", "16S RIBOSOMAL RNA", "16S RRNA", "RRNL", "LSU", "RRN16", "L-RRNA", "16S LARGE SUBUNIT RIBOSOMAL RNA", "LARGE SUBUNIT RIBOSOMAL RNA"],
         "18S": ["18S", "18S RIBOSOMAL RNA", "18S RRNA", "18S SMALL SUBUNIT RIBOSOMAL RNA", "SMALL SUBUNIT RIBOSOMAL RNA"],
         "28S": ["28S", "28S RIBOSOMAL RNA", "28S RRNA", "28S LARGE SUBUNIT RIBOSOMAL RNA", "LARGE SUBUNIT RIBOSOMAL RNA"],

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
stdnames =  ['12S', '16S', '18S', '28S', 'ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6', 'AK', 'CAD', 'EF1A', 'H3', 'RNApol', 'Wg']

misc = ["similar to cytochrome oxidase subunit I", "similar to cytochrome oxidase subunit 1", "similar to cytochrome c oxidase subunit I; COI", "similar to cytochrome c oxidase subunit I", "sequence contains partial cox1, tRNA-Leu and partial cox2 genes", "cox1", "coding region not determined; cytochrome oxidase subunit I", "similar to cytochrome oxidase subunit I; COI", "similar to cytochrome c oxidase", "similar to cytochrome oxidase subunit 1; 3' barcoding region; LCO-HCO"]

suborders = ['Adephaga', 'Polyphaga', 'Myxophaga', 'Archostemata']

unrec_genes = []
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
        output = genbank_metadata(rec)
        g = 0
        for feature in rec.features:
            type = feature.type
            if type not in ('CDS', 'rRNA'):
                continue
            else:
                name = get_feat_name(feature)
                stdname = ""
                for k, v in genes.items():
                    if name in v:
                        stdname = k
                        g += 1
                if stdname == '':
                    unrec_genes.append(name)
                    continue
                if args.mito:
                    if stdname not in mito:
                        continue
                if args.coi:
                    if stdname != 'COX1':
                        continue
                if stdname in cds:
                    if 'codon_start' in feature.qualifiers:
                        frame = feature.qualifiers["codon_start"][0]
                    else:
                        frame = ''
                        #print(f"Reading frame missing from record {rec.name}, {stdname}.")
                else:
                    frame = ''
            seq = ''
            seq = feature.extract(rec.seq)
            sequences.append(seq)
            output.update({"gene": stdname,
                           "length": len(seq),
                           "seq": seq,
                           "frame": frame})
            if output['txid'] in species:                              # If taxon ID in dict
                if stdname in species[output['txid']]:                 # If gene in dict for that taxon ID
                    species[output['txid']][stdname].append(output)    # Add gene info list to dict
                    x += 1
                else:
                    species[output['txid']][stdname] = [output]      # Otherwise add to dict with new key
                    x += 1
            else:
                species[output['txid']] = {stdname: [output]}      # Otherwise add to dict with new key
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

# Write CSV metadata file
gbids = []
with open("metadata.csv", "w") as file:
    writer = csv.writer(file)
    writer.writerow(
        ["ncbi_taxid", "genbank_accession", "bold_id", "bold_bin", "lab_id", "suborder", "infraorder", "superfamily", "family", 
        "subfamily", "tribe", "species", "country", "latitude", "longitude", "ref_authoer", "ref_title", "ref_journal"])
    for gene, records in longest.items():
        for rec in records:
            if rec['gbid'] not in gbids:
                gbids.append(rec['gbid'])
                row = [rec['txid'], '', '', '', rec['gbid']] + rec['taxonomy'] + [rec["spec"], rec['country'], rec['lat'], rec['long']] + rec['refs']
            writer.writerow(row)
    print("Metadata saved to metadata.csv")


# Write fastas
for gene, records in longest.items():
    file = open(f"{gene}.fasta", "w")
    x = 0
    y = 0
    for rec in records:
        if gene in rna:
            file.write(f">{rec['gbid']}\n{rec['seq']}\n")
            x += 1
        else:
            if rec['frame'] == '':
                if y == 0:
                    rf = open(f"{gene}.rf", "w")
                rf.write(f">{rec['gbid']}\n{rec['seq']}\n")
                y += 1
            else:
                file.write(f">{rec['gbid']};frame={rec['frame']}\n{rec['seq']}\n")
                x += 1

    print(f'{x} records written to {gene}.fasta')
    if y > 0:
        print(f'{y} records without reading frame written to {gene}.rf')


if args.accs:
    if len(nohits) > 0:
        print(f'\nNo requested genes found in the following records: {nohits}')

print("\nUnrecognised Genes")
counter = Counter(unrec_genes)
print(counter)
print('Misc Features')
print(misc_feature)
print("Other Feature Types")
print(other_type)