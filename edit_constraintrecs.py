
import csv

rowdict = {}
with open("constraint.csv") as file:
    csv_reader = csv.reader(file, delimiter=',')
    x = 0
    for row in csv_reader:
        if x == 0:
            x += 1
            continue
        else:
            x += 1
            tax = row[1]
            if tax in rowdict:
                rowdict[tax].append(row)
            else:
                rowdict[tax] = [row]

with open ("test.csv", "w") as file:
    writer = csv.writer(file)
    writer.writerow(["Accession", "Taxon ID", "Count", "Description", "Gene Count", "28S", "AK", "CAD", '12S', '16S', '18S', 'EF1A', 'H3',
         'Wg', 'ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6',
         "Date Late Modified", "Date Collected", "Domain", "Kingdom", "Superphylum", "Phylum",
         "Subphylum", "Class", "Subclass", "Infraclass", "Superorder", "Order", "Suborder", "Superfamily", "Family",
         "Subfamily", "Tribe", "Species", "Country", "Region", "Lat/Long", "Ref1 Author", "Ref1 Title", "Ref1 Journal",
         "Ref2 Author", "Ref2 Title", "Ref2 Journal", "Ref3 Author", "Ref3 Title", "Ref3 Journal"])
    for tax, rowlist in rowdict.items():
        n = 0
        m = 0
        for row in rowlist:
            if not all(g == "" for g in row[5:12]):
                n += 1
                #print(row[5:12])
            if not all(g == "" for g in row[13:27]):
                m += 1
        if m != 0 and n != 0:
            print(tax)
            for row in rowlist:
                writer.writerow(row)
