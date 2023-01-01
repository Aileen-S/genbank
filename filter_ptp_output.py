import csv

meta = {}
with open('metadata.csv') as file:
    metadata = csv.reader(file)
    for row in metadata:
        count = 0
        for r in row[3:20]:
            if r != '':
                count +=1
        meta[row[1]] = count


file = open('test.txt')
output = open('chosen.txt', 'w')

lines = file.readlines()
for line in lines:
    line = line.strip()
    line = line.split(',')
    if len(line) == 1:
        output.write(f'{line[0]}\n')
    else:
        txids = []
        chosen = line[0]
        for name in line:
            print(name)
            txid = name.split('_', 1)
            txids.append(txid[0])
            count = meta[txids[0]]
        for txid in txids:
            print(meta[txid])
            if meta[txid] > count:
                chosen = name
        print(f'Count = {count}')
        print(f'Chosen = {chosen}\n')
