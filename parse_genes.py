
from Bio import SeqIO
COX1_Records = []
for record in SeqIO.parse("test.gb", "gb"):
    for feature in record.features:
        for k, v in feature.qualifiers.items():     # Set up loop to search for key/value
            if k=="gene" and v==["COX1"]:           # Key = "gene", Value = gene name
                COX1_seq = record[feature.location.start:feature.location.end] # Extract sequence
                COX1_Records.append(COX1_seq)       # Add sequences to list

COX1full = []
for record in COX1_Records:
    if len(record.seq) >= 1300:     # From extracted COX1 records, refine to sequences over 1300bp
        COX1fasta = record.format("fasta")
        COX1full.append(COX1fasta)     # Add to new list
        print(COX1full, file = open("COX1fullseqs.txt", "w"))

