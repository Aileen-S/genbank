import argparse

parser = argparse.ArgumentParser(description="Create partition files from IQTREE output")
parser.add_argument("-i", "--input", type=str, help="in file")
parser.add_argument('-t', '--sequence_type', choices=['aa', 'nt'], help="Is input file from amino acid or nucleotide "
                                                                        "supermatrix.")

args = parser.parse_args()

file = open(args.input)
lines = file.readlines()
rrna = ['12S', '16S', '18S', '28S']

# AA gene partitions

if args.sequence_type == 'aa':
    outputaa = open(f'partitions_aa.txt', 'w')
    for line in lines:
        line = line.strip()
        split = line.split('=')
        start = split[0].split('/')
        if any(r in start[1] for r in rrna):
            start[0] = 'DNA'
        else:
            start[0] = 'mtART'
        outputaa.write(f'{start[0]}, {start[1]} = {split[1]}\n')

# NT Gene Partitions

if args.sequence_type == 'nt':
    outputnt = open(f'partitions_gene.txt', 'w')
    for line in lines:
        line = line.strip()
        split = line.split('=')
        start = split[0].split('/')
        start[0] = 'DNA'
        outputnt.write(f'{start[0]}, {start[1]} = {split[1]}\n')

    # NT Codon Partitions

    output123 = open(f'partitions_codon123.txt', 'w')
    output12 = open(f'partitions_codon12.txt', 'w')
    for line in lines:
        line = line.strip()
        split = line.split('=')
        start = split[0].split('/')
        start[0] = 'DNA'
        if any(r in start[1] for r in rrna):
            output123.write(f'{start[0]}, {start[1]} = {split[1]}\n')
            output12.write(f'{start[0]}, {start[1]} = {split[1]}\n')

        else:
            part = start[1].split('.')
            pos = split[1].split('-')
            output123.write(f'{start[0]}, {part[0]}.1 = {pos[0]}-{pos[1]}\\3\n')
            output12.write(f'{start[0]}, {part[0]}.1 = {pos[0]}-{pos[1]}\\3\n')
            output123.write(f'{start[0]}, {part[0]}.2 = {int(pos[0])+1}-{pos[1]}\\3\n')
            output12.write(f'{start[0]}, {part[0]}.2 = {int(pos[0])+1}-{pos[1]}\\3\n')
            output123.write(f'{start[0]}, {part[0]}.3 = {int(pos[0])+2}-{pos[1]}\\3\n')

