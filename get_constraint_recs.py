

# python3 get_constraint_recs.py -e mixedupvoyage@gmail.com -t Megadytes -n 1


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
parser = argparse.ArgumentParser(description="Search GenBank, retrieve gene sequences for records with at least specified number of genes from gene dict.", formatter_class=MultilineFormatter)
parser.add_argument("-t", "--taxon", type=str, help="Taxon of interest")
parser.add_argument("-n", "--number", type=int, help="Set minimum number of genes (from gene dict) required before record is saved")
#parser.add_argument("-g", "--gene", type=str, help="Gene(s) of interest. Format: gene1,gene2,gene3")
parser.add_argument("-e", "--email", type=str, help="Your email registered with NCBI")


# Start the actual script

args = parser.parse_args()         # Process input args from command line
#args = argparse.Namespace(taxon='Eretes', email='aileen.scott@nhm.ac.uk') # This is how I step through the script interactively
#Namespace(taxon='Eretes', mpc=True, email='aileen.scott@nhm.ac.uk', nuclear=False)


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




Entrez.email = args.email

handle = Entrez.egquery(term=args.taxon)
record = Entrez.read(handle)
for row in record["eGQueryResult"]:
    if row["DbName"]=="nuccore":
        count = int(row["Count"])
        print(f"{count} records found.")

species = {}
unrecgenes = set()
x = 0
y = 0
chunk = 1000
for start in range(0, count, chunk):
    print(f"Downloading GenBank records for taxon IDs {y} to {y+999}" if (y+999) < count else
      f"Downloading GenBank records for taxon IDs {y} to {count}")
    y += 1000
    handle = Entrez.esearch(db="nucleotide", term=args.taxon, retstart=start, retmax=chunk)
    record = Entrez.read(handle)
    gilist = record["IdList"]
    gi_str = ",".join(gilist)  # Join into string for efetch
    handle = Entrez.efetch(db="nucleotide", id=gi_str, rettype="gb", retmode="text")  # Get GenBanks
    record = SeqIO.parse(handle, "gb")

    for rec in record:
        z = 0
        gbid = rec.name
        #if args.taxon not in rec.annotations["taxonomy"]:
        #    continue

        # Get record output data
        db_xref = rec.features[0].qualifiers["db_xref"]
        for ref in db_xref:
            if "taxon" in ref:
                tax = "".join(filter(str.isdigit, ref))
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
            if "genes" in output:
                if stdname not in output["genes"]:
                    z += 1
                    output["genes"][stdname] = {"gene": stdname,
                                                "type": type,
                                                "length": len(seq),
                                                "seq": seq}
            else:
                z += 1
                output["genes"] = {stdname: {"gene": stdname,
                                             "type": type,
                                             "length": len(seq),
                                             "seq": seq}}
        if z >= args.number:
            output["gene count"] = z # Only save recs with at least specified number of genes
            #print(f"{gbid} has more than {args.number} gene(s)")
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


print(f"\n{str(x)} gene records saved to species dict")

print("\nUnrecognised Genes")
print(unrecgenes)

gen = ['18S', '28S', 'AK', 'CAD', 'EF1A', 'H3', 'RNApol', 'Wg', '12S', '16S', 'ATP6', 'ATP8',
       'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']

# Write CSV metadata file
with open("metadata.csv", "w") as file:     # Open output file
    writer = csv.writer(file)               # Name writer object
    writer.writerow(
        ["Accession", "Taxon ID", "Description", "Gene Count", '18S', "28S", "AK", "CAD", 'EF1A', 'H3', 'RNApol', 'Wg',
         '12S', '16S', 'ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6',
         "Date Late Modified", "Date Collected", "Domain", "Kingdom", "Superphylum", "Phylum",
         "Subphylum", "Class", "Subclass", "Infraclass", "Superorder", "Order", "Suborder", "Superfamily", "Family",
         "Subfamily", "Tribe", "Species", "Country", "Region", "Lat/Long", "Ref1 Author", "Ref1 Title", "Ref1 Journal",
         "Ref2 Author", "Ref2 Title", "Ref2 Journal", "Ref3 Author", "Ref3 Title", "Ref3 Journal"])

sequencedict = {}
file = open("metadata.csv", "a")
writer = csv.writer(file)
for rec in species.values():
    for gbid, output in rec.items():
        row = [output["gbid"], output["txid"], output["description"], output["gene count"]]
        for gene in gen:
            if gene in output["genes"]:
                row.append(output["genes"][gene]["length"])
            else:
                row.append("")
    # Also include gene count column
        row.append(output["rec date"])
        row.append(output["c date"])
        output["taxonomy"].extend([""] * (15 - len(output["taxonomy"])))
        if output["taxonomy"][14] == "Cybistrini":
            output["taxonomy"][13] = "Cybistrinae"
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



