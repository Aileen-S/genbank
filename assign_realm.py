from csv import writer
from csv import reader
import csv

file = open("biogeography.csv")
regions = {}

lines = file.readlines()
for line in lines:
    line = line.strip()
    realm, country = line.split(":")
    countries = country.split(",")
    regions[realm] = countries

#print(regions)

with open("metadata220614.csv", "r") as input, open("test.csv", "w") as output:
    reader = reader(input)
    writer = writer(output)
    next(reader, None)
    writer.writerow(["Accession", "Realm", "Country"])
    for row in reader:
        realm = ""
        for key in regions.keys():
            if row[21] in regions[key]:
                realm = key
        newrow = [row[0], realm, row[21]]
        writer.writerow(newrow)

# Make new csv instead with just accession and realm, then can copy and paste. Edit get_species list to include realm for next run.?


# Append a column in existing csv using csv.reader / csv.writer classes
def add_column_in_csv(input_file, output_file, transform_row):
    with open(input_file, 'r') as read_obj, \
            open(output_file, 'w', newline='') as write_obj:
        csv_reader = reader(read_obj)                   # Create a csv.reader object from the input file object
        csv_writer = writer(write_obj)                  # Create a csv.writer object from the output file object
        for row in csv_reader:                          # Read each row of the input csv file as list
            transform_row(row, csv_reader.line_num)     # Pass the list / row in the transform function to add column text for this row
            csv_writer.writerow(row)                    # Write the updated row / list to the output file
