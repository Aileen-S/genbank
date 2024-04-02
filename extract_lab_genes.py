import argparse
import csv
from Bio import Entrez
from Bio import SeqIO


# Function definitions

def get_feat_name(feat):
    featname = "unknown"
    nametags = ['gene', 'product', 'label', 'standard_name']  # Search these four keys for gene name
    if any(t in feat.qualifiers.keys() for t in nametags):
        for t in nametags:
            if t in feat.qualifiers.keys():
                featname = feat.qualifiers[t][0].upper()
                featname = featname.replace(' ', '')
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
#parser.add_argument('-m', '--mito', action='store_true', help='Save only mitochondrial protein-coding genes')
parser.add_argument('-i', '--fasta_id', choices=['gbid', 'txid', 'both'], help="Choose identifiers for output fastas. Default is gbid.")
parser.add_argument('-l', '--list', type=str, help="Limit to list of db_ids in file")


args = parser.parse_args()         # Process input args from command line

suborders = ['Adephaga', 'Polyphaga', 'Myxophaga', 'Archostemata']

mito = ['ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']

genes = {"12S": ["12S", "12S RIBOSOMAL RNA", "12S RRNA", "RRNS", "SSU"],
         "16S": ["16S", "16S RIBOSOMAL RNA", "16S RRNA", "RRNL", "LARGESUBUNITRIBOSOMALRNA", "LSU"],
         "ATP6": ['ATP SYNTHASE F0 SUBUNIT 6', 'APT6', 'ATP SYNTHASE A0 SUBUNIT 6', 'ATP SYNTHASE SUBUNIT 6', 'ATP SYNTHASE FO SUBUNIT 6', 'ATPASE6', 'ATPASE SUBUNIT 6', 'ATP6'],
         "ATP8": ['ATP SYNTHASE F0 SUBUNIT 8', 'APT8', 'ATP SYNTHASE A0 SUBUNIT 8', 'ATP SYNTHASE SUBUNIT 8', 'ATP SYNTHASE FO SUBUNIT 8', 'ATPASE8', 'ATPASE SUBUNIT 8', 'ATP8'],
         "COX1": ['CYTOCHROME C OXIDASE SUBUNIT 1', 'CYTOCHROME OXIDASE SUBUNIT I', 'CYTOCHROME C OXIDASE SUBUNIT I', 'COXI', 'CO1', 'COI', 'CYTOCHROME COXIDASE SUBUNIT I', 'CYTOCHROME OXIDASE SUBUNIT 1', 'CYTOCHROME OXYDASE SUBUNIT 1', 'COX1', 'CYTOCHROME OXIDASE I', 'CYTOCHROME OXIDASE C SUBUNIT I'],
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

ids = []
if args.list:
    file = open(args.list)
    lines = file.readlines()
    for line in lines:
        ids.append(line.strip())


# Search through GBIDs
species = {}
count = {}
with open(args.gb_file) as file:
    record = SeqIO.parse(file, "gb")
    for rec in record:
        if args.list:
            if rec.name not in ids:
                continue
        try:
            db_xref = rec.features[0].qualifiers["db_xref"]
            for ref in db_xref:
                if "taxon" in ref:  # Get NCBI taxon, rather than BOLD cross ref
                    txid = f'lab_{"".join(filter(str.isdigit, ref))}'  # Extract numbers from NCBI taxon value
        except KeyError:
            txid = ''
        try:
            spec = rec.annotations["organism"]
            # Replace the following characters: > < . ( ) ; : ' ,
            spec = spec.replace(">", "_").replace("<", "_").replace(".", "").replace('(', '_')\
                .replace(')', '_').replace(';', '_').replace(':', '_').replace("'", "").replace(',', '')
            specfasta = spec.replace(" ", "_")
        except KeyError:
            spec = ''
            specfasta = ''
        taxonomy = ['', '', '', '', '']
        try:
            for tax in rec.annotations["taxonomy"]:
                if tax in suborders: taxonomy[0] = tax
                if tax.endswith('oidea'): taxonomy[1] = tax
                if tax.endswith('idae'): taxonomy[2] = tax
                if tax.endswith('inae'): taxonomy[3] = tax
                if tax.endswith('ini'): taxonomy[4] = tax
            taxonomy.append(spec.split(' ')[0])
            fastatax = f"{taxonomy[2]}_{taxonomy[3]}_{taxonomy[4]}_{specfasta}"
        except KeyError:
            print(f'No taxonomy for {rec.name}')
        g = 0
        for feature in rec.features:
            type = feature.type
            if type in ('CDS', 'rRNA'):
                name = get_feat_name(feature)                       # Find gene name
                stdname = ""
                for k, v in genes.items():
                    if name in v:
                        stdname = k
                        g += 1
                    else:
                        continue
                #if stdname not in mito:
                #    continue
                if 'codon_start' in feature.qualifiers:
                    frame = feature.qualifiers["codon_start"]
                else:
                    frame = ''
                seq = feature.extract(rec.seq)
                output = {"gene": stdname,
                          "gbid": rec.name,
                          "txid": txid,
                          "description": rec.description,
                          "spec": spec,
                          "taxonomy": taxonomy,
                          "fastatax": fastatax,
                          "length": len(seq),
                          "seq": seq,
                          "frame": frame}
                if stdname in species:
                    species[stdname].append(output)
                else:
                    species[stdname] = [output]
        count[rec.name] = g

meta = open("metadata.csv", "w")
writer = csv.writer(meta)
writer.writerow(["db_id", "count"])
for k, v in count.items():
    writer.writerow([k, v])

rna = ['12S', '16S', '18S', '28S']



for gene, records in species.items():
    if gene in rna:
        file = open(f"{gene}.fasta", "w")
        x = 0
        for rec in records:
            fasta_id = f">{rec['gbid']};frame={rec['frame'][0]}\n{rec['seq']}\n"
            file.write(fasta_id)
            x += 1
            print(f'{x} records written to {gene}.fasta')

    else:
        file = open(f"{gene}.fasta", "w")
        rf = open(f"{gene}.fa", "w")
        x = 0
        y = 0
        for rec in records:
            if rec['frame'] == '':
                fasta_id = f">{rec['gbid']}\n{rec['seq']}\n"
            else:
                fasta_id = f">{rec['gbid']};frame={rec['frame'][0]}\n{rec['seq']}\n"

            if 'frame' in fasta_id:
                file.write(fasta_id)
                x += 1
            else:
                rf.write(fasta_id)
                y += 1
        print(f'{x} records written to {gene}.fasta')
        print(f'{y} records without reading frame written to {gene}.fa')
