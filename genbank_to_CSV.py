# Parse GenBank file to get CSV file with accession numbers and taxonomy
# "Domain", "Kingdom", "Superphylum", "Phylum", "Subphylum", "Class", "Subclass", "Infraclass", "Superorder", "Order", "Suborder", "Superfamily", "Family", "Subfamily", "Tribe"
# Need to extract db_xref for csv file

import csv
from Bio import SeqIO

with open("COX1_longest.csv", "w") as file:             # Open output file
    csvfile = csv.writer(file)                          # Name writer object
    csvfile.writerow(["Accession", "taxon_id", "Species", "Domain", "Kingdom", "Superphylum", "Phylum", "Subphylum", "Class", "Subclass", "Infraclass", "Superorder", "Order", "Suborder", "Superfamily", "Family", "Subfamily", "Tribe", "1", "2", "3"])  # Specify column names
    for record in SeqIO.parse("COX1_longest.gb", "gb"): # Open GenBank file
        Acc = record.id[0:-2]                                 # Get accessions
        for feature in record.features:
            for k, v in feature.qualifiers.items():
                if k == "db_xref":                      # Search for db_xref (taxon ID)
                    str_id = str(v)                     # Save as string
                    if "taxon" in str_id:
                        xref = str_id[8:-2]             # Extract numbers. Original format is eg ['taxon:207455']
        Spe = record.annotations["organism"]            # Get genus/species
        Tax = record.annotations["taxonomy"]            # Get higher taxonomy list
        Row = [Acc, xref, Spe]                          # New list with accession and species names
        for t in Tax:
            Row.append(t)                               # Add desired taxon levels to list
        csvfile.writerow(Row)                           # Write accession/taxonomy list to row

