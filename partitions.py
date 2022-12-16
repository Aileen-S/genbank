# python3 partitions.py -i test.txt -t aa -o nexus

import argparse

parser = argparse.ArgumentParser(description="Create partition files from IQTREE output")
parser.add_argument("-i", "--input", type=str, help="in file")
parser.add_argument('-t', '--sequence_type', choices=['aa', 'nt'], help="Is input file from amino acid or nucleotide "
                                                                        "supermatrix.")
parser.add_argument('-o', '--output', choices=['raxml', 'nexus'], help="Choose output file type")


args = parser.parse_args()

file = open(args.input)
lines = file.readlines()
rrna = ['12S', '16S', '18S', '28S']


# AA gene partitions

if args.sequence_type == 'aa':
    outputaa = open(f'partitions_aa.txt', 'w')
    if args.output == 'nexus':
        outputaa.write('\#nexus\n'
                       'begin sets;\n')
    for line in lines:
        line = line.strip()
        split = line.split('=')         # split[0] eg all_aa/12S.fasta, split[1] eg 1-436
        start = split[0].split('/')     # start[0] eg all_aa, start[1] eg 12s.fasta
        if any(r in start[1] for r in rrna):    # Check for rRNA genes
            start[0] = 'DNA'
        else:
            start[0] = 'mtART'
        if args.output == 'raxml':
            outputaa.write(f'{start[0]}, {start[1]}={split[1]}\n')
        if args.output == 'nexus':
            outputaa.write(f'  charset {start[1]}={split[1]}\n')

    if args.output == 'nexus':
        outputaa.write('end;')


# NT Gene Partitions

if args.sequence_type == 'nt':
    outputnt = open(f'partitions_gene.txt', 'w')
    if args.output == 'nexus':
        outputnt.write('\#nexus\n'
                       'begin sets;\n')
    for line in lines:
        line = line.strip()
        split = line.split('=')
        start = split[0].split('/')
        start[0] = 'DNA'
        if args.output == 'raxml':
            outputnt.write(f'{start[0]}, {start[1]}={split[1]}\n')
        if args.output == 'nexus':
            outputnt.write(f'  charset {start[1]}={split[1]}\n')


    # NT Codon Partitions

    output123 = open(f'partitions_codon123.txt', 'w')
    output12 = open(f'partitions_codon12.txt', 'w')
    if args.output == 'nexus':
        output123.write('\#nexus\n'
                       'begin sets;\n')
        output12.write('\#nexus\n'
                       'begin sets;\n')
    for line in lines:
        line = line.strip()
        split = line.split('=')
        start = split[0].split('/')
        start[0] = 'DNA'
        if any(r in start[1] for r in rrna):
            if args.output == 'raxml':
                output123.write(f'{start[0]}, {start[1]}={split[1]}\n')
                output12.write(f'{start[0]}, {start[1]}={split[1]}\n')
            if args.output == 'nexus':
                output123.write(f'  charset {start[1]}={split[1]}\n')
                output12.write(f'  charset {start[1]}={split[1]}\n')

        else:
            part = start[1].split('.')
            pos = split[1].split('-')
            if args.output == 'raxml':
                output123.write(f'{start[0]}, {part[0]}.1 ={pos[0]}-{pos[1]}\\3\n')
                output12.write(f'{start[0]}, {part[0]}.1 ={pos[0]}-{pos[1]}\\3\n')
                output123.write(f'{start[0]}, {part[0]}.2 = {int(pos[0])+1}-{pos[1]}\\3\n')
                output12.write(f'{start[0]}, {part[0]}.2 = {int(pos[0])+1}-{pos[1]}\\3\n')
                output123.write(f'{start[0]}, {part[0]}.3 = {int(pos[0])+2}-{pos[1]}\\3\n')
            if args.output == 'nexus':
                output123.write(f'  charset {part[0]}.1 ={split[1]}\n')
                output123.write(f'  charset {part[0]}.2 ={split[1]}\n')
                output123.write(f'  charset {part[0]}.3 ={split[1]}\n')
                output12.write(f'  charset {part[0]}.1 ={split[1]}\n')
                output12.write(f'  charset {part[0]}.2 ={split[1]}\n')
    if args.output == 'nexus':
        outputnt.write('end;')
        output12.write('end;')
        output123.write('end;')




