from Bio import SeqIO
import urllib

def loadnamevariants():
    output = {}
    url = "https://raw.githubusercontent.com/tjcreedy/constants/master/gene_name_variants.txt"
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
fh = open("COX1fullseqa.fasta", "w")        # Write new file
n = 0                                       # Start count
COX1full = []                               # Make empty list
for record in COX1_Records:                 # Process records found from genbank search
    if len(record.seq) >= 1000:             # Keep only sequences over set length
        n += 1                              # Count kept sequences
        COX1fasta = record.format("fasta")  # Format as fasta
        COX1full.append(COX1fasta)          # Add to new list
        print(COX1full, file = open("COX1fullseqs.txt", "w"))   # Print to file

print(str(n) + " records saved")