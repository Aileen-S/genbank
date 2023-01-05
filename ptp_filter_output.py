import csv

meta = {}
with open('metadata.csv') as file:  # output from ptp_get_metadata.py
    metadata = csv.reader(file)
    for row in metadata:
        count = 0
        for r in row[3:20]:
            if r != '':
                count +=1
        meta[row[1]] = count

lab = ['50515', '50516', '107801', '107841', '107843', '107861', '107887']
found = []

file = open('output.txt')     # output from [grep "^ " bPTP_221219.PTPhSupportPartition.txt > output.txt]
chosen = open('chosen.txt', 'w')

names = []
lines = file.readlines()
for line in lines:
    line = line.strip()
    line = line.split(',')
    if len(line) == 1:
        chosen.write(f'{line[0]}\n')
    else:
        txids = []
        best = line[0]
        for name in line:
            #print(name)
            if '_sp_' in name or '_BM' in name:
                print('sp or bm')
                print(name)
                txid = name.split('_', 1)
                if txid[0] in lab:
                    found.append(txid[0])
                else:
                    txids.append(txid[0])
                    count = meta[txids[0]]
            else:
                chosen.write(f'{name}\n')

        for txid in txids:
            #print(meta[txid])
            if meta[txid] > count:
                names.append(name)

for name in names:
    chosen.write(f'{name}\n')

        #print(f'Count = {count}')
        #print(f'Chosen = {best}\n')

print(f'Lab records: {found}')