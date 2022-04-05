from Bio import Entrez
import pandas as pd
from Bio import SeqIO

Entrez.email = "aileen.scott@nhm.ac.uk"

# Search to see how many records there are
handle = Entrez.egquery(term="Dytiscidae AND NAD4")                             
record = Entrez.read(handle)
for row in record["eGQueryResult"]:
    if row["DbName"]=="nuccore":    # Don't understand when I should use nuccore/nucleotide
        print(row["Count"])

# Get the list of accession numbers.
# Default is 20 records, changed retmax to maximum allowed.
handle = Entrez.esearch(db="nuccore", term="Dytiscidae AND NAD4", retmax=100000)
record = Entrez.read(handle)
accessions = record["IdList"]
print(record["Count"])

# Use accession numbers to download GenBank records
# Did this because couldn't get efetch to read my list of accessions in the next step.
# Couldn't work out how to use retstart to return >10,000 (for COX1)
acc_str = ",".join(accessions)
handle = Entrez.efetch(db="nucleotide", id=acc_str, rettype="gb", retmode="text")
print(handle.read(), file=open("GBNAD4.gb", "w"))

# Check you have the right number of sequences
GBfile = open("GBNAD4.gb")
GBfile = GBfile.read()
count = GBfile.count("LOCUS")
print(str(count) + " records found")




# Find longest sequences

with open("GBNAD4.gb") as file: 
    Accs = []                                       # Make empty lists
    Species = []
    Lengths = []
    for rec in SeqIO.parse(file, "genbank"):        # Start for loop
        Accs.append(rec.id)                         # Add accessions to Accs list
        Species.append(rec.annotations["organism"]) # Add species name to Species list
        Lengths.append(len(rec.seq))                # Add sequence length to Lengths list


# Combine lists
zipped = list(zip(Accs, Species, Lengths))


# Make dataframe, naming columns
df = pd.DataFrame(zipped, columns = ["Accessions", "Species", "Sequence Length"])

print("Initial Dataframe")
print("Dataframe has " + str(len(df)) + " rows.")
print(df.head(5))


# Sort alphabetically by species.
# inpace=True overwrites existing dataframe. Default is false.
df.sort_values(by=['Species'], inplace=True)


# Sort by species (a-z), then sequence length (descending order)
df.sort_values(by = ["Species", "Sequence Length"], ascending = [True, False], inplace=True)

print("\nSorted Dataframe")
print(df.head(5))


# Delete duplicate species, keeping first (longest sequence)
df.drop_duplicates(         # Remove duplicates
    subset = ["Species"],   # Specify column
    keep = "first",         # Keep first occurrance
    inplace=True)           # Overwrite dataframe

print("\nDulpicates removed, longest sequence for each species retained.")
print("Dataframe has " + str(len(df)) + " rows.")
print(df.head(5))


# Extract column with accession numbers
Acc=(df["Accessions"].tolist())

# Write it to file with no quotes or commas.
with open("/users/aileenscott/Desktop/GenBank/Barcodes/TestNAD4.txt", mode="w") as file:
    file.write("\n".join(Acc))

