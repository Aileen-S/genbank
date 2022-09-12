import csv

Ambiguous = ['Mexico', 'Angola', 'Mozambique', 'Namibia', 'Algeria', 'Botswana', 'Chad', 'Democratic Republic of the Congo', 'Egypt', 'Liyba', 'Mauritania', 'Niger', 'Nigeria', 'Mali', 'Afghanistan', 'China', 'Saudi Arabia', 'Russia', 'Canada', 'USA', 'Mexico', 'Papua New Guinea', 'Argentina', 'Bolivia', 'Brazil', 'Chile', 'Peru']

file = open("countries.txt")
regions = {}

lines = file.readlines()
for line in lines:
    line = line.strip()
    region, other = line.split(";")
    subregion, countries = other.split(":")
    countries = countries.split(",")
    if region in regions:
        regions[region][subregion] = countries
    else:
        regions[region] = {subregion: countries}


x = 0
with open("metadata_all.csv") as input, open("test.csv", "w") as output:
    reader = csv.reader(input)
    for row in reader:
        x += 1
        if x == 1:
            writer = csv.writer(output)
            writer.writerow(row)
        else:
            for region, sub in regions.items():
                for subregion, countries in sub.items():
                    if row[26] in countries:
                        row[24] = region
                        row[25] = subregion
                if row[26] in Ambiguous:
                    row[23] = 'area data needed'
            writer.writerow(row)