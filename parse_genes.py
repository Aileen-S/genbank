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

x = 0
y = 0
COX1_Records = []
for record in SeqIO.parse("COX1_longest.gb", "gb"):     # Start with genabnk file
    for feature in record.features:                     # Examine each record feature
        for k, v in feature.qualifiers.items():         # Search key:values pairs
            if v in genenames:          # Not sure of the right syntax for this bit
                y += 1
                stdname = genenames[v]
                COX1_seq = record[feature.location.start:feature.location.end] # Extract sequence
                if len(COX1_seq) >= 100:                    # Check sequence length
                    x += 1
                    COX1_fasta = COX1_seq.format("fasta")   # Save as fasta
                    COX1_Records.append(COX1_fasta)         # Keep sequences over set length
                print(COX1_Records, file=open("COX1fullseqs.fasta", "w"))  # Print to file

print(str(y) + " records found")
print(str(x) + " records saved")