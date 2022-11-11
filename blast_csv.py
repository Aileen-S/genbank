from Bio.Blast import NCBIXML
import argparse
import csv

parser = argparse.ArgumentParser(description='BLAST search input fasta. Save results list and interim XML file.')
parser.add_argument("-i", "--input", type=str, help="Input xml.")
parser.add_argument("-o", "--output", type=str, help="Output csv file")
args = parser.parse_args()


print(f'Writing output to {args.output}')
output = open(f"{args.output}", "w")
writer = csv.writer(output)
writer.writerow(['Query', 'Hit', 'e value', 'Length', 'Start'])

range = {1:30}
for n in range:
    result_handle = open(f"1vas{n}.fa.xml")
    blast_records = NCBIXML.parse(result_handle)

    E_VALUE_THRESH = 0.04
    x = 0
    for blast_record in blast_records:
        x += 1
        y = 0
        print(f'\nBLAST SEARCH {x}: {blast_record.query}\n')
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if y == 5:
                    break
                if hsp.expect < E_VALUE_THRESH:
                    y +=1
                    #print(hsp.)
                    print(f"hit: {alignment.title}")
                    writer.writerow([blast_record.query, alignment.title, hsp.expect, alignment.length, hsp.sbjct_start])
                    print(f"length:   {alignment.length}")
                    print(f"e value:  {hsp.expect}")
                    print(f"starts:   {hsp.sbjct_start}")
                    print(hsp.query[0:99] + "...")
                    print(hsp.match[0:99] + "...")
                    print(hsp.sbjct[0:99] + "...")
            if y == 5:
                break
