import argparse
from Bio import Entrez
from Bio import SeqIO
import csv

parser = argparse.ArgumentParser(description="Rename sequences in fasta file from CSV.")
parser.add_argument("-i", "--input", type=str, help="Input BLAST result ( -outfmt '6 staxids sacc sseq')")
parser.add_argument("-o", "--output", type=str, help="Output csv with old and new names")
parser.add_argument("-m", "--metadata", type=str, help="Optional output CSV file to get metadata")
parser.add_argument("-e", "--email", type=str, help="Your email address for NCBI")

args = parser.parse_args()

records = {}
x = 0
with open(args.input, "r") as file:
    data = file.readlines()
    for line in data:
        x += 1
        line = line.strip()
        record = line.split("\t") # Record format 'TXID \t GBID \t seq'
        txid = record[0]
        gbid = record[1]
        seq = record[2]
        if txid not in records:
            records[txid] = {gbid: seq}
        else:
            records[txid][gbid] = seq
print(f'{x} records in {args.input}')
print(f'{len(records)} unique taxon IDs')

longest = {}
gbids = []
for txid, recs in records.items():
    length = 0
    for gbid, seq in recs.items():
        if len(seq) > length:
            longest[txid] = {gbid: seq}
            chosen = gbid
    longest[txid] = {gbid: seq}
    gbids.append(chosen)

with open(args.output, 'w') as output:
    for txid, rec in longest.items():
        for gbid, seq in rec.items():
            output.write(f'>{txid}_{gbid}\n{seq}\n')
print(f'Saved longest sequences to {args.output}')
print(f'Searching NCBI for metadata')

if args.metadata:
    Entrez.email = args.email
    suborders = ['Adephaga', 'Polyphaga', 'Myxophaga', 'Archostemata']
    metadata = []
    accstr = ",".join(gbids)
    handle = Entrez.efetch(db="nucleotide", id=accstr, rettype="gb", retmode="text")
    record = SeqIO.parse(handle, "gb")
    for rec in record:
        db_xref = rec.features[0].qualifiers["db_xref"]
        bold = ''
        for ref in db_xref:
            if "taxon" in ref:  # Get NCBI taxon, rather than BOLD cross ref
                txid = "".join(filter(str.isdigit, ref))  # Extract numbers from NCBI taxon value
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
            continue
        output = {"gbid": rec.name,
                  "txid": txid,
                  "description": rec.description,
                  "spec": rec.annotations["organism"],
                  "rec date": rec.annotations["date"],
                  "c date": c_date,
                  "taxonomy": taxonomy,
                  "fastatax": fastatax,
                  "country": country,
                  "region": region,
                  "latlon": latlon,
                  "lat": lat,
                  "long": long,
                  "refs": refs}
        metadata.append(output)

    with open(args.metadata, "w") as output:
        writer = csv.writer(output)  # Name writer object
        writer.writerow(
            ["FastaID", "Taxon ID", "Accession", "Species", "Suborder", "Superfamily", "Family", "Subfamily", "Tribe", 'Genus', "Description", "Date Late Modified",
             "Date Collected", "Country", "Region", "Lat/Long", "Lat", "Long", "Ref1 Author", "Ref1 Title",
             "Ref1 Journal", "Ref2 Author",
             "Ref2 Title", "Ref2 Journal", "Ref3 Author", "Ref3 Title", "Ref3 Journal"])
        for output in metadata:
            row = [f"{output['txid']}_{output['fastatax']}", output["txid"], output["gbid"],
                   output["spec"]]
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
print(f'Saved metadata to {args.metadata}')