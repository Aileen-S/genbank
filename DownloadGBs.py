from Bio import Entrez
Entrez.email = "mixedupvoyage@gmail.com"



handle = Entrez.egquery(term="Agabus")
record = Entrez.read(handle)
for row in record["eGQueryResult"]:
    if row["DbName"]=="nuccore":
        count = int(row["Count"])
print(count)

chunk = 1000
for start in range(0, count, chunk):
    handle = Entrez.esearch(db="nucleotide", term="Agabus", retstart=start, retmax=chunk)
    record = Entrez.read(handle)
    gi_list = record["IdList"]
    #print(record["Count"])
    print(gi_list)

# Use GIs to download GenBank records
gi_str = ",".join(gi_list)
for start in range(0, count, chunk):
    handle = Entrez.efetch(db="nucleotide", id=gi_str, rettype="gb", retmode="text", retstart=start, retmax=chunk) # Need to figure out retstart
    print(handle.read(), file=open("Agabus.gb", "w"))

