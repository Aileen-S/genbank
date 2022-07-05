# Run interactively
# From /mbl/share/workspaces/groups/voglerlab/MMGdatabase/gbmaster_2022-06-27

# Use file with list of dbids to combine individual genbank records into one file, then convert to fasta
file = open("/home/ailes/mitogenomes/mito_dbids_220704.txt")
output = open("/home/ailes/mitogenomes/mito_genbanks_220704.gb", "a")

lines = file.readlines()
for line in lines:
    line = line.strip()
    name = line + ".gb"
    record = open(name)
    record = record.read()
    output.write(record)

from Bio import SeqIO
genbank = SeqIO.parse("/home/ailes/mitogenomes/mito_genbanks_220704.gb", "gb")
fasta = open("/home/ailes/mitogenomes/mito_genbanks_220704.fa", "a")

for record in genbank:
    fasta.write(f">{record.name}\n{record.seq}\n")


# From ~/mitogenomes

# Get records from genbank search with >5 genes and append to fasta file
from Bio import Entrez
from Bio import SeqIO
Entrez.email = ""

file = open("genbank_5andover.txt")
lines = file.readlines()
for line in lines:
    line = line.strip()
    handle = Entrez.efetch(db="nucleotide", id=line, rettype="gb", retmode="text")
    record = SeqIO.parse(handle, "gb")
    for rec in record:
        fasta = open("mito_220704.fa", "a")
        fasta.write(f">{rec.name}\n{rec.seq}\n")

