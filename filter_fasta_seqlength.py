# Find sequences over certain length in fasta file and write to new fasta file

from Bio import SeqIO

COX1000 = []
for record in SeqIO.parse("8_gbnt_raw/COX1.fa", "fasta"):
    if len(record.seq) >= 1000:     # From extracted COX1 records, refine to sequences over set length
        COX1000.append(record)     # Add to new list
        SeqIO.write(COX1000, "example.fasta", "fasta")