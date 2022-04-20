# Parse GenBank file to get CSV file with accession numbers and taxonomy
# "Domain", "Kingdom", "Superphylum", "Phylum", "Subphylum", "Class", "Subclass", "Infraclass", "Superorder", "Order",
# "Suborder", "Superfamily", "Family", "Subfamily", "Tribe"
# Need to extract db_xref for csv file

import csv
from Bio import SeqIO

with open("COX1_longest.csv", "w") as file:             # Open output file
    csvfile = csv.writer(file)                          # Name writer object
    csvfile.writerow(["Accession", "Species", "Domain", "Kingdom", "Superphylum", "Phylum", "Subphylum", "Class", "Subclass", "Infraclass", "Superorder", "Order", "Suborder", "Superfamily", "Family", "Subfamily", "Tribe", "1", "2", "3"])  # Specify column names
    for record in SeqIO.parse("COX1_longest.gb", "gb"): # Open GenBank file
        Acc = record.name                      # Get accessions (record.name gives accession, record.id gives version)
        Spe = record.annotations["organism"]            # Get genus/species
        Tax = record.annotations["taxonomy"]            # Get higher taxonomy list
        Row = [Acc, Spe]                                # New list with accession and species names
        for t in Tax:
            Row.append(t)                               # Add desired taxon levels to list
        csvfile.writerow(Row)                           # Write accession/taxonomy list to row

print(record.name)