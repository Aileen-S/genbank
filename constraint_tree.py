import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Filter fasta file from list of IDs")
parser.add_argument("-i", "--input", type=str, help="Input fasta")
parser.add_argument("-o", "--output", type=str, help="Output constraint tree")

args = parser.parse_args()

Dytiscidae_families = {'Coptotominae': [],
                       'Lancetinae': [],
                       'Agabinae': [],
                       'Colymbetinae': [],
                       'Cybistrini': [],
                       'Cybistrini2': [],
                       'Dytiscinae': [],
                       'Dytiscinae2': [],
                       'Copelatinae': [],
                       'Matinae': [],
                       'Laccophilinae': [],
                       'Hydroporinae': []
                       }

Dytiscidae = '((Coptotominae,Lancetinae),(Laccophilinae,(((Agabinae,Colymbetinae),(Cybistrini2,Dytiscinae2)),(Copelatinae,(Matinae,Hydroporinae)))));'

records = SeqIO.parse(args.input, "fasta")
x = 0
for rec in records:
    for key in Dytiscidae_families.keys():
        if key in rec.id:
            Dytiscidae_families[key].append(rec.id)

for taxon in list(Dytiscidae_families['Dytiscinae']):
    if 'Cybist' not in taxon:
        Dytiscidae_families['Dytiscinae2'].append(taxon)

Dytiscidae_families.pop('Dytiscinae')
Dytiscidae_families['Cybistrini2'] = Dytiscidae_families.pop('Cybistrini')

for key, value in Dytiscidae_families.items():
    if value == []:
        print(f'Warning, no taxa found for {key}')
    else:
        print(f'{len(value)} taxa for {key}')

for key, value in Dytiscidae_families.items():
    taxa = ','.join(value)
    Dytiscidae = Dytiscidae.replace(key, f'({taxa})')

with open(args.output, 'w') as file:
    file.write(Dytiscidae)

