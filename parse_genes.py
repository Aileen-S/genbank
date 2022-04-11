
from Bio import SeqIO
import urllib

def loadnamevariants():
    output = {}
    url = "https://raw.githubusercontent.com/tjcreedy/genenames/master/gene_name_variants.txt"
    for line in urllib.request.urlopen(url):
        line = line.decode('utf-8').strip()
        name = line.split(";")[0]
        annotype = line.split(":")[0].split(";")[1]
        variants = line.split(":")[1].split(",")
        for v in variants:
            for g in ['', ' ']:
                v = v.replace(g, '')
                for s in ['',' GENE', ' '+annotype.upper()]:
                    output[v+s] = name
    return(output)

# Get the gene name variants
genenames = loadnamevariants()

COX1_Records = []
for record in SeqIO.parse("COX1_longest.gb", "gb"):
    for feature in record.features:
        if feature.id in genenames:
            stdname = genenames[feature.id]
            if stdname==["COX1"]:
                COX1_seq = record[feature.location.start:feature.location.end] # Extract sequence
                COX1_Records.append(COX1_seq)       # Add sequences to list

# Filter out longest sequences
COX1full = []
for record in COX1_Records:
    if len(record.seq) >= 600:             # From extracted COX1 records, refine to sequences over 1300bp
        COX1fasta = record.format("fasta")
        COX1full.append(COX1fasta)          # Add to new list
        print(COX1full, file = open("COX1fullseqs.txt", "w"))

# Check how many sequences
GBfile = open("COX1fullseqs.txt")
GBfile = GBfile.read()
count = GBfile.count(">")
print(str(count) + " records saved")