# Run interactively
# From /mbl/share/workspaces/groups/voglerlab/MMGdatabase/gbmaster_2022-06-27

# Use file with list of dbids to combine individual genbank records into one file.
file = open("/home/ailes/mitogenomes/mito_dbids_220704.txt")
output = open("/home/ailes/mitogenomes/mito_genbanks_220704.gb", "a")

lines = file.readlines()
for line in lines:
    line = line.strip()
    name = line + ".gb"
    record = open(name)
    record = record.read()
    output.write(record)


# From ~/mitogenomes

# Get records from genbank search with >5 genes and append to genbank file
from Bio import Entrez
Entrez.email = ""

file = open("genbank_5andover.txt")
lines = file.readlines()
for line in lines:
    line = line.strip()
    handle = Entrez.efetch(db="nucleotide", id=line, rettype="gb", retmode="text")
    print(handle.read(), file=open("mito_genbanks_220704.gb", "a"))
