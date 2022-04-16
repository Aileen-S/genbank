# Parse GenBank file to get CSV file with accession numbers and taxonomy

import csv
from Bio import SeqIO

with open("COX1_longest.CSV", "w") as file:             # Open output file
    csvfile = csv.writer(file)                          # Name writer object
    csvfile.writerow(["accession", "species", "family", "subfamily", "tribe"])  # Specify column names
    for record in SeqIO.parse("COX1_longest.gb", "gb"): # Open GenBank file
        Acc = record.id                                 # Get accessions
        Spe = record.annotations["organism"]
        Tax = record.annotations["taxonomy"]            # Get taxonomy list
        Row = [Acc, Spe]                                # New list with accession in first column
        n = 0                                           # Start counter
        for t in Tax:
            n += 1                                      # Count loops
            if n >= 13 and n <= 15:                     # Skip taxon levels we don't need
                Row.append(t)                           # Add desired taxon levels to list
        csvfile.writerow(Row)                           # Write accession/taxonomy list to row