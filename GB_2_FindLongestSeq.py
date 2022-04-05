# Get longest sequence for each species

# Remember to import panda and SeqIO
import pandas as pd
from Bio import SeqIO


with open("COX1_genbank_download.gb") as file:      # Specify location and name file
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
print(df.head(10))


# Sort alphabetically by species.
# inpace=True overwrites existing dataframe. Default is false.
df.sort_values(by=['Species'], inplace=True)


# Sort by species (a-z), then sequence length (descending order)
df.sort_values(by = ["Species", "Sequence Length"], ascending = [True, False], inplace=True)

print("\nSorted Dataframe")
print(df.head(10))


# Delete duplicate species, keeping first (longest sequence)
df.drop_duplicates(         # Remove duplicates
    subset = ["Species"],   # Specify column
    keep = "first",         # Keep first occurrance
    inplace=True)           # Overwrite dataframe

print("\nDulpicates removed, longest sequence retained.")
print("Dataframe has " + str(len(df)) + " rows.")
print(df.head(10))


# Extract column with accession numbers
Acc=(df["Accessions"].tolist())
print(Acc)

# Write to file
print(Acc, file=open("Accessions_List.txt", mode="w"))

# Write it to file with no quotes or commas.
with open("Accessions_Simple.txt", mode="w") as file:
    file.write("\n".join(Acc))

