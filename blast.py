#python3 blast.py -i test.txt -o blast.txt

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import argparse

parser = argparse.ArgumentParser(description='BLAST search input fasta. Save results list and interim XML file.')
parser.add_argument("-i", "--input", type=str, help="Input fasta")
parser.add_argument("-o", "--output", type=str, help="Output txt file")
args = parser.parse_args()


fasta_string = open(args.input).read()
print(f"Running BLAST for sequence(s) in {args.input}")
result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string)
print(f"Writing results to {args.input}.xml")
with open(f"{args.input}.xml", "w") as out_handle:
    out_handle.write(result_handle.read())

print(f"Parsing {args.input}.xml")
result_handle = open(f"{args.input}.xml")
blast_records = NCBIXML.parse(result_handle)

print(f'Writing output to {args.output}')
output = open(args.output, "w")
E_VALUE_THRESH = 0.04
x = 0
for blast_record in blast_records:
    x += 1
    y = 0
    output.write(f'\nBLAST SEARCH {x}: {blast_record.query}\n\n')
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if y == 10:
                break
            if hsp.expect < E_VALUE_THRESH:
                y +=1
                output.write(f"Hit {y}\n")
                #print(hsp.)
                output.write(f"sequence: {alignment.title}\n")
                output.write(f"length:   {alignment.length}\n")
                output.write(f"e value:  {hsp.expect}\n")
                output.write(f"starts:   {hsp.sbjct_start}\n")
                output.write(hsp.query[0:99] + "...\n")
                output.write(hsp.match[0:99] + "...\n")
                output.write(hsp.sbjct[0:99] + "...\n")
        if y == 10:
            break

    if y == 0:
        output.write(f'\n{blast_record.query}\n'
              f'No hits with e value < 0.04.\n\n')
        print("{blast_record.query}\n"
              "No hits with e value < 0.04.")
    else:
        print(f'{blast_record.query}\n'
              f'Saved {y} records with e value < 0.04')

