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
                       'Dytiscinae': [],
                       'Copelatinae': [],
                       'Matinae': [],
                       'Laccophilinae': [],
                       'Hydroporinae': []}

Dytiscidae = '((Coptotominae,Lancetinae),(Laccophilinae,(((Agabinae,Colymbetinae),(Cybistrini,Dytiscinae)),(Copelatinae,(Matinae,Hydroporinae)))));'

records = SeqIO.parse(args.input, "fasta")
x = 0
for rec in records:
    for key in Dytiscidae_families.keys():
        if key in rec.id:
            id = rec.id.replace('.', '').replace('-', '')
            Dytiscidae_families[key].append(rec.id)

for key, value in Dytiscidae_families.items():
    if value == []:
        print(f'Warning, no taxa found for {key}')

for key, value in Dytiscidae_families.items():
    taxa = ','.join(value)
    Dytiscidae = Dytiscidae.replace(key, f'({taxa})')
    #print(taxa)

file = open(args.output, 'w')
file.write(Dytiscidae)

print(Dytiscidae_families['Cybistrini'])