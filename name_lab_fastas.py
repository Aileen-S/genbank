import argparse

parser = argparse.ArgumentParser(description="Get list of accessions or taxon IDs from fasta or list of sequence names.")
parser.add_argument("-i", "--input", type=str, help="in file")
args = parser.parse_args()


label = \
    {'MIZA00211': '50515_labMIZA00211_Dytiscidae_Beetle_1096',
     'MIZA00183': '50516_labMIZA00183_Dytiscidae_Dytiscinae_Cybistrini_Cybister_24',
     'MIZA00004': '107801_labMIZA00004_Dytiscidae_Dytiscinae_Eretini_Eretes_83',
     'BIOD00670': '107841_labBIOD00670_Dytiscidae_Copelatinae_Copelatus_SF_Dyt1',
     'BIOD00322': '107843_labBIOD00322_Dytiscidae_Hydroporinae_Bidessini_MS965',
     'MIZA00010': '107861_labMIZA00010_Dytiscidae_Dytiscinae_Hydaticini_Hydaticus_14',
     'MIZA00023': '107887_labMIZA00023_Dytiscidae_Hydroporinae_Hyphydrini_Hyphydrus_3',
     'MIZA00288': '107891_labMIZA00288_Laccophilinae_Laccophilini_Laccophilus_43',
     'GBDL01546': '254283_labGBDL01546_Dytiscidae_Hydroporinae_Hydroporini_Hydroporus_palustris',
     'GBDL02028': '940473_labGBDL02028_Dytiscidae_Dytiscinae_Dytiscini_Dytiscus_sharpi',
     'GBDL01697': '1205547_labGBDL01697_Dytiscidae_Copelatinae_Copelatus_sp._COP02',
     'GBDL02687': '2811616_labGBDL02687_Dytiscidae_Hydroporinae_Hydroporini_Hydroporus_dobrogeanus_complex_sp._IBE<ESP>_AV49',
     'MIZA00232': '50515_labMIZA00232_Dytiscidae_Beetle_1096'}

refs = ['MIZA00211', 'MIZA00183', 'MIZA00004', 'BIOD00670', 'BIOD00322', 'MIZA00010', 'MIZA00023', 'MIZA00288',
        'GBDL01546', 'GBDL02028', 'GBDL01697', 'GBDL02687', 'MIZA00232']

file = open(args.input)
lines = file.readlines()
output = open(args.input, 'w')
for line in lines:
    for k, v in label.items():
        if k in line:
            line = line.replace(k, v)
    output.write(line)
