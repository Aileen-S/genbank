import argparse
from Bio import Entrez
from Bio import SeqIO
import csv
import time


parser = argparse.ArgumentParser(description="Rename sequences in fasta file from CSV.")
parser.add_argument("-i", "--input", type=str, help="Input BLAST result ( -outfmt '6 staxids sacc sseq')")
parser.add_argument("-o", "--output", type=str, help="Output fasta")
parser.add_argument("-m", "--metadata", type=str, help="Optional output CSV file to get metadata")
parser.add_argument("-e", "--email", type=str, help="Your email address for NCBI, to get metadata")

args = parser.parse_args()


# Search GenBank, wait and try again if 'too many requests' error
def search_genbank(ids):
    for attempt in range(1, 10):
        try:
            handle = Entrez.efetch(db="nucleotide", id=ids, rettype="gb", retmode="text")
            results = SeqIO.parse(handle, "gb")
            return results
        except Entrez.HTTPError:
            print("HTTP error fetching records; retry in 20 seconds")
            time.sleep(20)
    print(f"Failed to retrieve records after 10 attempts.")
    return None

Entrez.email = args.email

# Save records in dict goruped by TXID
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

        # Check for multiple TXIDs
        if ';' in txid:
            results = search_genbank(gbid)
            for r in results:
                db_xref = r.features[0].qualifiers["db_xref"]
                for ref in db_xref:
                    if "taxon" in ref:  # Get NCBI taxon, rather than BOLD cross ref
                        txid = "".join(filter(str.isdigit, ref))  # Extract numbers from NCBI taxon value
        if txid not in records:
            records[txid] = {gbid: seq}
        else:
            records[txid][gbid] = seq
print(f'{x} records in {args.input}')
print(f'{len(records)} unique NCBI taxon IDs')

# Find longest sequence for each TXID
print(f'Finding longest sequence for each NCBI taxon ID')
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

# Write fasta
with open(args.output, 'w') as output:
    for txid, rec in longest.items():
        for gbid, seq in rec.items():
            output.write(f'>{txid}_{gbid}\n{seq}\n')
print(f'Saved longest sequences to {args.output}')

# Get metadata
if args.metadata:
    print(f'Searching NCBI for metadata')
    suborders = ['Adephaga', 'Polyphaga', 'Myxophaga', 'Archostemata']
    metadata = []
    accstr = ",".join(gbids)
    results = search_genbank(accstr)
    for rec in results:
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
            first = rec.annotations['references'][0]
            refs.append(first.authors)
            refs.append(first.title)
            refs.append(first.journal)
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
        writer.writerow(["FastaID", "Taxon ID", "Accession", "Species", "Suborder", "Superfamily", "Family",
                         "Subfamily", "Tribe", 'Genus', "Description", "Date Late Modified", "Date Collected",
                         "Country", "Region", "Lat/Long", "Lat", "Long", "Ref1 Author", "Ref1 Title", "Ref1 Journal"])
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