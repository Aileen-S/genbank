from Bio import SeqIO
import argparse
import csv

parser = argparse.ArgumentParser(description="Rename sequences in fasta file from CSV.")
parser.add_argument("-i", "--input", type=str, help="Input fasta")
parser.add_argument("-n", "--names", type=str, help="List of new names (with old name as part of name")
parser.add_argument("-c", "--csv", type=str, help="CSV file, new names in first column, old names in second")
parser.add_argument("-o", "--output", type=str, help="Output fasta file")
parser.add_argument("-r", "--renamed", type=str, help="Output csv with old and new names")
args = parser.parse_args()

recs = {}

if args.names:
    new = []
    with open(args.names) as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip()
            new.append(line)

    records = SeqIO.parse(args.input, "fasta")
    for rec in records:
        new_id = rec.id
        if ';frame==' in rec.id:
            r_id, frame = rec.id.split(';')
            for n in new:
                if r_id in n:
                    new_id = f'{n};{frame}'
                elif rec.id == n:
                    new_id = f'{n};{frame}'

        else:
            for n in new:
                if rec.id in n:
                    new_id = n
                elif rec.id == n:
                    new_id = n
        if new_id in recs:
            recs[new_id][rec.id] = rec.seq
        else:
            recs[new_id] = {rec.id: rec.seq}

if args.csv:
    # new names in first column, old names in second
    meta = {}
    with open(args.csv) as file:
        metadata = csv.reader(file)
        for row in metadata:
            if row [0] != '' and row[0] != 'NA':
                meta[row[1]] = row[0]

    records = SeqIO.parse(args.input, "fasta")
    output = open(args.output, 'w')
    for rec in records:
        new_id = rec.id
        if ';frame=' in rec.id:
            r_id, frame = rec.id.split(';')
            for k, v in meta.items():
                if r_id == k:
                    new_id = f'{v};{frame}'
        else:
            for k, v in meta.items():
                if rec.id == k:
                    new_id = v
        if new_id in recs:
            recs[new_id][rec.id] = rec.seq
        else:
            recs[new_id] = {rec.id: rec.seq}


# Check for duplicate names
# Save longest sequence for duplicates
selected = {}
removed = {}

for new, old in recs.items():
    max_len = 0
    max_rec = ''
    for k, v in old.items():
        if len(v) > max_len:
            max_rec = {k: v}
            max_len = len(v)
        else:
            removed[new] = k
    selected[new] = max_rec

#print(f'Removed shorter sequences from duplicate IDs: {removed}')

if args.renamed:
    with open(args.renamed, 'w') as id_list:
        writer = csv.writer(id_list)
        writer.writerow(['Old ID', 'New ID', 'Removed'])
        for new, rec in selected.items():
            for old, seq in rec.items():
                writer.writerow([old, new])
        for new, old in removed.items():
            writer.writerow([old, new, 'yes'])
    print(f'Saved original IDs to {args.renamed}')


output = open(args.output, 'w')
for new, rec in selected.items():
    for old, seq in rec.items():
        output.write(f'>{new}\n{seq}\n')

print(f'Saved renamed fasta to {args.output}')
