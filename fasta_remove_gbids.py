import argparse

parser = argparse.ArgumentParser(description="Remove accessions from fasta.")
parser.add_argument("-i", "--input", type=str, help="in file")
parser.add_argument("-o", "--output", type=str, help="out file")


args = parser.parse_args()

file = open(args.input)
output = open(args.output, 'w')
lines = file.readlines()
for line in lines:
    try:
        if line.startswith('>'):
            name = line.split("_")
            if "NC" in name:
                del name[1:3]
            else:
                if any(n.isdigit() for n in name[1]):
                    del name[1]
            line = '_'.join(name)
    except IndexError:
        pass
    output.write(line)


