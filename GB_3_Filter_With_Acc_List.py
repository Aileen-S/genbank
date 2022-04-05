from Bio import Entrez, SeqIO

Entrez.email = "aileen.scott@nhm.ac.uk"

# Found this in the biopython manual and only edited slightly.
# Don't completely understand how it works.

# Specify files to read/write
input_file = "GBNAD4.gb"        # Genbank file with all sequences
id_file = "NAD4.txt"     # List of accessions for selected sequences
output_file = "NAD.gb"         # Output genbank file for selected sequences

with open(id_file) as id_handle:
    wanted = set(line.rstrip("\n").split(None, 1)[0] for line in id_handle)
    
print("Found %i unique identifiers in %s" % (len(wanted), id_file))

records = (r for r in SeqIO.parse(input_file, "gb") if r.id in wanted)
count = SeqIO.write(records, output_file, "gb")

print("Saved %i records from %s to %s" % (count, input_file, output_file))

if count < len(wanted):
    print("Warning %i IDs not found in %s" % (len(wanted) - count, input_file))
