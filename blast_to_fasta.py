import argparse

parser = argparse.ArgumentParser(description="Rename sequences in fasta file from CSV.")
parser.add_argument("-i", "--input", type=str, help="Input BLAST result ( -outfmt '6 staxids sacc sseq')")
parser.add_argument("-o", "--output", type=str, help="Output csv with old and new names")
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
        if gbid not in records:
            records[txid] = {gbid: seq}
        else:
            records[txid][gbid] = seq
print(f'{x} records in {args.input}')
print(f'{len(records)} unique taxon IDs')

longest = {}
for txid, recs in records.items():
    length = 0
    for gbid, seq in recs.items():
        if len(seq) > length:
            chosen = gbid
    longest[txid] = {gbid: seq}

with open(args.output, 'w') as output:
    for txid, rec in longest.items():
        for gbid, seq in rec.items():
            output.write(f'>{txid}_{gbid}\n{seq}\n')
print(f'Saved longest sequences to {args.output}')