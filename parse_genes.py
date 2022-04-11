from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature

from Bio import SeqIO
COX1_Records = []
for record in SeqIO.parse("test.gb", "gb"):
    for feature in record.features:
        for k, v in feature.qualifiers.items():
            if k=="gene" and v==["COX1"]:
                if 
                COX1_seq = record[feature.location.start:feature.location.end]
                COX1_Records.append(COX1_seq)

for record in COX1_Records:
    print(feature.qualifiers.items)