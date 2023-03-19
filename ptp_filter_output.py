import csv

meta = {}
with open('metadataptp.csv') as file:  # output from ptp_get_metadata.py
    metadata = csv.reader(file)
    for row in metadata:
        count = 0
        for r in row[3:20]:
            if r != '':
                count += 1
        meta[row[1]] = count    # Save count of genes for each taxon ID

#lab = ['50515', '50516', '107801', '107841', '107843', '107861', '107887']
#found = []

file = open('output.txt')     # output from [grep "^ " bPTP_221219.PTPhSupportPartition.txt > output.txt]
chosen = open('chosen.txt', 'w')

names = {}
x = 0
most = []
lines = file.readlines()
for line in lines:
    line = line.strip()
    line = line.split(',')
    if len(line) == 1:
        chosen.write(f'{line[0]}\n') # Write single taxon lines to file
    else:
        x += 1
        txids = []
        for name in line:
            txid = name.split('_', 1)[0] # Get taxon IDs for multi taxon lines
            txids.append(txid)
            names[txid] = name          # Save fasta IDs
            count = meta[txids[0]]      # Save number of genes from first record on line
            top = txids[0]              # Save taxon ID from first record on line
        for txid in txids:
            if meta[txid] > count:      # Replace count and taxon ID if another record has more genes
                count = meta[txid]
                top = txid
        most.append(top)    # Save taxon ID with most genes

for txid in most:
    chosen.write(f'{names[txid]}\n')    # Write fasta ID of longest

#for name in names:
 #   chosen.write(f'{name}\n')

        #print(f'Count = {count}')
        #print(f'Chosen = {best}\n')

#print(f'Lab records: {found}')