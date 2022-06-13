from Bio import Entrez
Entrez.email = "mixedupvoyage@gmail.com"



handle = Entrez.egquery(term="Dytiscidae")
record = Entrez.read(handle)
for row in record["eGQueryResult"]:
    if row["DbName"]=="nuccore":
        count = int(row["Count"])
print(count)

chunk = 1000
gi_list = []
for start in range(0, count, chunk):
    handle = Entrez.esearch(db="nucleotide", term="Dytiscidae", retstart=start, retmax=chunk)
    record = Entrez.read(handle)
    for gi in record["IdList"]:
        gi_list.append(gi)
    print(len(gi_list))

# Use GIs to download GenBank records
gi_str = ",".join(gi_list)
for start in range(0, count, chunk):
    handle = Entrez.efetch(db="nucleotide", id=gi_str, rettype="gb", retmode="text", retstart=start, retmax=chunk) # Need to figure out retstart
    print(handle.read(), file=open("Dytiscidae.gb", "a"))

GBfile = open("Dytiscidae.gb")
GBfile = GBfile.read()
count = GBfile.count("LOCUS")
print(str(count) + " records found")