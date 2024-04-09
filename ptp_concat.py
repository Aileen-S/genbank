
import csv
import argparse

# Argument parser
parser = argparse.ArgumentParser(description="Process PTP output file.")
parser.add_argument("-i", "--input", type=str, help="mPTP output txt file")
parser.add_argument("-o", "--output", type=str, help="output taxon list with PTP species numbers")
args = parser.parse_args()


# Save mPTP delimited species lists
ptp_meta = {}
x = 0
y = 0
with open(args.input) as file:
    lines = file.readlines()
    for line in lines:
        y += 1
        if 'Number of delimited species' in line:
            print(line)
        x += 1
        if y < 8:
            continue
        if 'Species ' in line:
            spec_no = line.strip().replace(':', '')
            spec_no = spec_no.split(' ')[1]
            spec_no = str(spec_no).zfill(4)
            spec_no = f'T{spec_no}'
            continue
        line = line.strip()
        if line != '':
            if spec_no in ptp_meta:
                ptp_meta[spec_no].append(line)
            else:
                ptp_meta[spec_no] = [line]

# Write PTP metadata
with open(args.output, 'w') as file:
    writer = csv.writer(file)
    writer.writerow(['PTP Species', 'Taxon'])
    print(f'Writing output to {args.output}')
    for k, values in ptp_meta.items():
        print(k, values)
        x += 1
        for v in values:
            writer.writerow([k, v])
