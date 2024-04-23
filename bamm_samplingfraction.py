import argparse

subgen = {'Agabus': ['Acatodes', 'Gaurodytes'],
          'Platynectes': ['Agametrus', 'Australonectes', 'Gueorguievtes', 'Leuronectes'],
          'Cybister': ['Megadytoides', 'Melanectes', 'Neocybister'],
          'Megadytes': ['Bifurcitus', 'Paramegadytes', 'Trifurcitus'],
          'Acilus': ['Homoeolytrus'],
          'Hydaticus': ['Prodaticus'],
          'Clypeodytes': ['Hypoclypeus', 'Paraclypeus'],
          'Clemnius': ['Cyclopius'],
          'Hygrotus': ['Coelambus', 'Heroceras', 'Herophydrus', 'Hyphoporus', 'Leptolambus'],
          'Rhantus': ['Anisomera', 'Senilites'],
          'Aglymbus': ['Rugosus'],
          'Exocelina': ['Papuadytes'],
          'Paroster': ['Terradessus']}

sub = {'Agabinae': {'total': 426, 'count': 0, 'fraction': 0},
       'Colymbetinae': {'total': 142, 'count': 0, 'fraction': 0},
       'Copelatinae': {'total': 791, 'count': 0, 'fraction': 0},
       'Coptotominae': {'total': 5, 'count': 0, 'fraction': 0},
       'Cybistrinae': {'total': 129, 'count': 0, 'fraction': 0},
       'Dytiscinae': {'total': 254, 'count': 0, 'fraction': 0},
       'Hydrodytinae': {'total': 4, 'count': 0, 'fraction': 0},
       'Hydroporinae': {'total': 2373, 'count': 0, 'fraction': 0},
       'Laccophilinae': {'total': 460, 'count': 0, 'fraction': 0},
       'Lancetinae': {'total': 22, 'count': 0, 'fraction': 0},
       'Matinae': {'total': 9, 'count': 0, 'fraction': 0}}

tribe = {'Aciliini': {'total': 70, 'count': 0, 'fraction': 0},
         'Agabetini': {'total': 2, 'count': 0, 'fraction': 0},
         'Agabini': {'total': 338, 'count': 0, 'fraction': 0},
         'Aubehydrini': {'total': 2, 'count': 0, 'fraction': 0},
         'Bidessini': {'total': 752, 'count': 0, 'fraction': 0},
         'Colymbetini': {'total': 142, 'count': 0, 'fraction': 0},
         'Copelatini': {'total': 791, 'count': 0, 'fraction': 0},
         'Coptotomini': {'total': 5, 'count': 0, 'fraction': 0},
         'Cybistrini': {'total': 129, 'count': 0, 'fraction': 0},
         'Dytiscini': {'total': 29, 'count': 0, 'fraction': 0},
         'Eretini': {'total': 4, 'count': 0, 'fraction': 0},
         'Hydaticini': {'total': 149, 'count': 0, 'fraction': 0},
         'Hydrodytini': {'total': 4, 'count': 0, 'fraction': 0},
         'Hydrotrupini': {'total': 2, 'count': 0, 'fraction': 0},
         'Hydrovatini': {'total': 218, 'count': 0, 'fraction': 0},
         'Hygrotini': {'total': 140, 'count': 0, 'fraction': 0},
         'Hyphydrini': {'total': 395, 'count': 0, 'fraction': 0},
         'Laccophilini': {'total': 458, 'count': 0, 'fraction': 0},
         'Laccornellini': {'total': 40, 'count': 0, 'fraction': 0},
         'Laccornini': {'total': 10, 'count': 0, 'fraction': 0},
         'Lancetini': {'total': 22, 'count': 0, 'fraction': 0},
         'Matini': {'total': 9, 'count': 0, 'fraction': 0},
         'Methlini': {'total': 42, 'count': 0, 'fraction': 0},
         'Pachydrini': {'total': 14, 'count': 0, 'fraction': 0},
         'Platynectini': {'total': 86, 'count': 0, 'fraction': 0},
         'Vatellini': {'total': 59, 'count': 0, 'fraction': 0},
         'Deronectina': {'total': 197, 'count': 0, 'fraction': 0},
         'Hydroporina': {'total': 282, 'count': 0, 'fraction': 0},
         'Siettitiina': {'total': 59, 'count': 0, 'fraction': 0},
         'Sternopriscina': {'total': 155, 'count': 0, 'fraction': 0}}

gen = {'Agabinus': {'total': 2, 'count': 0, 'fraction': 0}, 'Agabus': {'total': 177, 'count': 0, 'fraction': 0},
       'Hydronebrius': {'total': 4, 'count': 0, 'fraction': 0}, 'Ilybiosoma': {'total': 17, 'count': 0, 'fraction': 0},
       'Ilybius': {'total': 71, 'count': 0, 'fraction': 0}, 'Platambus': {'total': 67, 'count': 0, 'fraction': 0},
       'Hydrotrupes': {'total': 2, 'count': 0, 'fraction': 0}, 'Andonectes': {'total': 2, 'count': 0, 'fraction': 0},
       'Platynectes': {'total': 83, 'count': 0, 'fraction': 0}, 'Bunites': {'total': 1, 'count': 0, 'fraction': 0},
       'Caperhantus': {'total': 1, 'count': 0, 'fraction': 0}, 'Carabdytes': {'total': 10, 'count': 0, 'fraction': 0},
       'Colymbetes': {'total': 22, 'count': 0, 'fraction': 0}, 'Hoperius': {'total': 1, 'count': 0, 'fraction': 0},
       'Meladema': {'total': 4, 'count': 0, 'fraction': 0}, 'Melanodytes': {'total': 1, 'count': 0, 'fraction': 0},
       'Meridiorhantus': {'total': 5, 'count': 0, 'fraction': 0}, 'Nartus': {'total': 2, 'count': 0, 'fraction': 0},
       'Neoscutopterus': {'total': 2, 'count': 0, 'fraction': 0}, 'Rhantus': {'total': 93, 'count': 0, 'fraction': 0},
       'Agaporomorphus': {'total': 12, 'count': 0, 'fraction': 0}, 'Aglymbus': {'total': 12, 'count': 0, 'fraction': 0},
       'Capelatus': {'total': 1, 'count': 0, 'fraction': 0}, 'Copelatus': {'total': 460, 'count': 0, 'fraction': 0},
       'Exocelina': {'total': 209, 'count': 0, 'fraction': 0}, 'Lacconectus': {'total': 80, 'count': 0, 'fraction': 0},
       'Liopterus': {'total': 2, 'count': 0, 'fraction': 0}, 'Madaglymbus': {'total': 15, 'count': 0, 'fraction': 0},
       'Coptotomus': {'total': 5, 'count': 0, 'fraction': 0}, 'Austrodytes': {'total': 2, 'count': 0, 'fraction': 0},
       'Cybister': {'total': 97, 'count': 0, 'fraction': 0}, 'Megadytes': {'total': 21, 'count': 0, 'fraction': 0},
       'Onychohydrus': {'total': 2, 'count': 0, 'fraction': 0}, 'Regimbartina': {'total': 1, 'count': 0, 'fraction': 0},
       'Spencerhydrus': {'total': 2, 'count': 0, 'fraction': 0}, 'Sternhydrus': {'total': 4, 'count': 0, 'fraction': 0},
       'Acilius': {'total': 13, 'count': 0, 'fraction': 0}, 'Aethionectes': {'total': 6, 'count': 0, 'fraction': 0},
       'Graphoderus': {'total': 12, 'count': 0, 'fraction': 0}, 'Rhantaticus': {'total': 1, 'count': 0, 'fraction': 0},
       'Sandracottus': {'total': 17, 'count': 0, 'fraction': 0},
       'Thermonectus': {'total': 20, 'count': 0, 'fraction': 0},
       'Tikoloshanes': {'total': 1, 'count': 0, 'fraction': 0}, 'Notaticus': {'total': 2, 'count': 0, 'fraction': 0},
       'Dytiscus': {'total': 27, 'count': 0, 'fraction': 0}, 'Hyderodes': {'total': 2, 'count': 0, 'fraction': 0},
       'Eretes': {'total': 4, 'count': 0, 'fraction': 0}, 'Hydaticus': {'total': 149, 'count': 0, 'fraction': 0},
       'Hydrodytes': {'total': 3, 'count': 0, 'fraction': 0},
       'Microhydrodytes': {'total': 1, 'count': 0, 'fraction': 0},
       'Africodytes': {'total': 5, 'count': 0, 'fraction': 0}, 'Allodessus': {'total': 5, 'count': 0, 'fraction': 0},
       'Amarodytes': {'total': 11, 'count': 0, 'fraction': 0}, 'Anodocheilus': {'total': 22, 'count': 0, 'fraction': 0},
       'Belladessus': {'total': 5, 'count': 0, 'fraction': 0}, 'Bidessodes': {'total': 20, 'count': 0, 'fraction': 0},
       'Bidessonotus': {'total': 36, 'count': 0, 'fraction': 0}, 'Bidessus': {'total': 52, 'count': 0, 'fraction': 0},
       'Borneodessus': {'total': 1, 'count': 0, 'fraction': 0}, 'Brachyvatus': {'total': 4, 'count': 0, 'fraction': 0},
       'Clypeodytes': {'total': 39, 'count': 0, 'fraction': 0}, 'Comaldessus': {'total': 1, 'count': 0, 'fraction': 0},
       'Crinodessus': {'total': 1, 'count': 0, 'fraction': 0}, 'Fontidessus': {'total': 7, 'count': 0, 'fraction': 0},
       'Geodessus': {'total': 2, 'count': 0, 'fraction': 0}, 'Gibbidessus': {'total': 8, 'count': 0, 'fraction': 0},
       'Glareadessus': {'total': 2, 'count': 0, 'fraction': 0}, 'Hemibidessus': {'total': 6, 'count': 0, 'fraction': 0},
       'Huxelhydrus': {'total': 1, 'count': 0, 'fraction': 0}, 'Hydrodessus': {'total': 33, 'count': 0, 'fraction': 0},
       'Hydroglyphus': {'total': 92, 'count': 0, 'fraction': 0}, 'Hypodessus': {'total': 6, 'count': 0, 'fraction': 0},
       'Incomptodessus': {'total': 1, 'count': 0, 'fraction': 0},
       'Kakadudessus': {'total': 1, 'count': 0, 'fraction': 0}, 'Leiodytes': {'total': 31, 'count': 0, 'fraction': 0},
       'Limbodessus': {'total': 79, 'count': 0, 'fraction': 0}, 'Liodessus': {'total': 48, 'count': 0, 'fraction': 0},
       'Microdessus': {'total': 1, 'count': 0, 'fraction': 0},
       'Neobidessodes': {'total': 10, 'count': 0, 'fraction': 0},
       'Neobidessus': {'total': 29, 'count': 0, 'fraction': 0},
       'Neoclypeodytes': {'total': 29, 'count': 0, 'fraction': 0},
       'Novadessus': {'total': 1, 'count': 0, 'fraction': 0}, 'Pachynectes': {'total': 3, 'count': 0, 'fraction': 0},
       'Papuadessus': {'total': 2, 'count': 0, 'fraction': 0}, 'Peschetius': {'total': 12, 'count': 0, 'fraction': 0},
       'Petrodessus': {'total': 1, 'count': 0, 'fraction': 0}, 'Platydytes': {'total': 4, 'count': 0, 'fraction': 0},
       'Pseuduvarus': {'total': 2, 'count': 0, 'fraction': 0}, 'Rompindessus': {'total': 1, 'count': 0, 'fraction': 0},
       'Sharphydrus': {'total': 4, 'count': 0, 'fraction': 0}, 'Sinodytes': {'total': 1, 'count': 0, 'fraction': 0},
       'Spanglerodessus': {'total': 1, 'count': 0, 'fraction': 0},
       'Tepuidessus': {'total': 2, 'count': 0, 'fraction': 0},
       'Trogloguignotus': {'total': 1, 'count': 0, 'fraction': 0},
       'Tyndallhydrus': {'total': 1, 'count': 0, 'fraction': 0}, 'Uvarus': {'total': 66, 'count': 0, 'fraction': 0},
       'Yola': {'total': 48, 'count': 0, 'fraction': 0}, 'Yolina': {'total': 12, 'count': 0, 'fraction': 0},
       'Zimpherus': {'total': 1, 'count': 0, 'fraction': 0}, 'Amurodytes': {'total': 1, 'count': 0, 'fraction': 0},
       'Boreonectes': {'total': 11, 'count': 0, 'fraction': 0}, 'Clarkhydrus': {'total': 10, 'count': 0, 'fraction': 0},
       'Deronectes': {'total': 61, 'count': 0, 'fraction': 0}, 'Deuteronectes': {'total': 2, 'count': 0, 'fraction': 0},
       'Hornectes': {'total': 1, 'count': 0, 'fraction': 0}, 'Iberonectes': {'total': 1, 'count': 0, 'fraction': 0},
       'Larsonectes': {'total': 1, 'count': 0, 'fraction': 0}, 'Leconectes': {'total': 1, 'count': 0, 'fraction': 0},
       'Mystonectes': {'total': 5, 'count': 0, 'fraction': 0}, 'Nebrioporus': {'total': 57, 'count': 0, 'fraction': 0},
       'Nectoboreus': {'total': 3, 'count': 0, 'fraction': 0}, 'Nectomimus': {'total': 1, 'count': 0, 'fraction': 0},
       'Nectoporus': {'total': 9, 'count': 0, 'fraction': 0}, 'Neonectes': {'total': 3, 'count': 0, 'fraction': 0},
       'Oreodytes': {'total': 14, 'count': 0, 'fraction': 0}, 'Scarodytes': {'total': 11, 'count': 0, 'fraction': 0},
       'Stictotarsus': {'total': 3, 'count': 0, 'fraction': 0}, 'Trichonectes': {'total': 1, 'count': 0, 'fraction': 0},
       'Zaitzevhydrus': {'total': 1, 'count': 0, 'fraction': 0}, 'Haideoporus': {'total': 1, 'count': 0, 'fraction': 0},
       'Heterosternuta': {'total': 14, 'count': 0, 'fraction': 0},
       'Hydrocolus': {'total': 12, 'count': 0, 'fraction': 0}, 'Hydroporus': {'total': 191, 'count': 0, 'fraction': 0},
       'Neoporus': {'total': 39, 'count': 0, 'fraction': 0},
       'Sanfilippodytes': {'total': 25, 'count': 0, 'fraction': 0},
       'Ereboporus': {'total': 1, 'count': 0, 'fraction': 0}, 'Etruscodytes': {'total': 1, 'count': 0, 'fraction': 0},
       'Graptodytes': {'total': 23, 'count': 0, 'fraction': 0}, 'Iberoporus': {'total': 4, 'count': 0, 'fraction': 0},
       'Lioporeus': {'total': 2, 'count': 0, 'fraction': 0}, 'Metaporus': {'total': 2, 'count': 0, 'fraction': 0},
       'Porhydrus': {'total': 4, 'count': 0, 'fraction': 0}, 'Psychopomporus': {'total': 1, 'count': 0, 'fraction': 0},
       'Rhithrodytes': {'total': 6, 'count': 0, 'fraction': 0}, 'Siettitia': {'total': 2, 'count': 0, 'fraction': 0},
       'Stictonectes': {'total': 12, 'count': 0, 'fraction': 0}, 'Stygoporus': {'total': 1, 'count': 0, 'fraction': 0},
       'Antiporus': {'total': 15, 'count': 0, 'fraction': 0}, 'Barretthydrus': {'total': 3, 'count': 0, 'fraction': 0},
       'Brancuporus': {'total': 2, 'count': 0, 'fraction': 0}, 'Carabhydrus': {'total': 10, 'count': 0, 'fraction': 0},
       'Chostonectes': {'total': 6, 'count': 0, 'fraction': 0}, 'Megaporus': {'total': 11, 'count': 0, 'fraction': 0},
       'Necterosoma': {'total': 12, 'count': 0, 'fraction': 0}, 'Paroster': {'total': 52, 'count': 0, 'fraction': 0},
       'Sekaliporus': {'total': 2, 'count': 0, 'fraction': 0},
       'Sternopriscus': {'total': 29, 'count': 0, 'fraction': 0}, 'Tiporus': {'total': 13, 'count': 0, 'fraction': 0},
       'Laodytes': {'total': 1, 'count': 0, 'fraction': 0}, 'Siamoporus': {'total': 1, 'count': 0, 'fraction': 0},
       'Tassilodytes': {'total': 1, 'count': 0, 'fraction': 0}, 'Hydrovatus': {'total': 215, 'count': 0, 'fraction': 0},
       'Queda': {'total': 3, 'count': 0, 'fraction': 0}, 'Clemnius': {'total': 8, 'count': 0, 'fraction': 0},
       'Hygrotus': {'total': 132, 'count': 0, 'fraction': 0}, 'Agnoshydrus': {'total': 8, 'count': 0, 'fraction': 0},
       'Allopachria': {'total': 47, 'count': 0, 'fraction': 0}, 'Andex': {'total': 1, 'count': 0, 'fraction': 0},
       'Anginopachria': {'total': 3, 'count': 0, 'fraction': 0}, 'Coelhydrus': {'total': 1, 'count': 0, 'fraction': 0},
       'Darwinhydrus': {'total': 1, 'count': 0, 'fraction': 0},
       'Desmopachria': {'total': 136, 'count': 0, 'fraction': 0},
       'Dimitshydrus': {'total': 1, 'count': 0, 'fraction': 0}, 'Hovahydrus': {'total': 4, 'count': 0, 'fraction': 0},
       'Hydropeplus': {'total': 2, 'count': 0, 'fraction': 0}, 'Hyphovatus': {'total': 3, 'count': 0, 'fraction': 0},
       'Hyphydrus': {'total': 140, 'count': 0, 'fraction': 0}, 'Microdytes': {'total': 47, 'count': 0, 'fraction': 0},
       'Primospes': {'total': 1, 'count': 0, 'fraction': 0}, 'Canthyporus': {'total': 38, 'count': 0, 'fraction': 0},
       'Laccornellus': {'total': 2, 'count': 0, 'fraction': 0}, 'Laccornis': {'total': 10, 'count': 0, 'fraction': 0},
       'Celina': {'total': 34, 'count': 0, 'fraction': 0}, 'Methles': {'total': 8, 'count': 0, 'fraction': 0},
       'Heterhydrus': {'total': 5, 'count': 0, 'fraction': 0}, 'Pachydrus': {'total': 9, 'count': 0, 'fraction': 0},
       'Derovatellus': {'total': 42, 'count': 0, 'fraction': 0}, 'Vatellus': {'total': 17, 'count': 0, 'fraction': 0},
       'Kuschelydrus': {'total': 1, 'count': 0, 'fraction': 0}, 'Morimotoa': {'total': 3, 'count': 0, 'fraction': 0},
       'Phreatodessus': {'total': 2, 'count': 0, 'fraction': 0},
       'Typhlodessus': {'total': 1, 'count': 0, 'fraction': 0}, 'Agabetes': {'total': 2, 'count': 0, 'fraction': 0},
       'Africophilus': {'total': 19, 'count': 0, 'fraction': 0},
       'Australphilus': {'total': 2, 'count': 0, 'fraction': 0}, 'Hamadiana': {'total': 1, 'count': 0, 'fraction': 0},
       'Japanolaccophilus': {'total': 1, 'count': 0, 'fraction': 0},
       'Laccodytes': {'total': 12, 'count': 0, 'fraction': 0}, 'Laccomimus': {'total': 13, 'count': 0, 'fraction': 0},
       'Laccophilus': {'total': 289, 'count': 0, 'fraction': 0}, 'Laccoporus': {'total': 1, 'count': 0, 'fraction': 0},
       'Laccosternus': {'total': 2, 'count': 0, 'fraction': 0}, 'Napodytes': {'total': 1, 'count': 0, 'fraction': 0},
       'Neptosternus': {'total': 98, 'count': 0, 'fraction': 0},
       'Philaccolilus': {'total': 13, 'count': 0, 'fraction': 0},
       'Philaccolus': {'total': 5, 'count': 0, 'fraction': 0}, 'Philodytes': {'total': 1, 'count': 0, 'fraction': 0},
       'Lancetes': {'total': 22, 'count': 0, 'fraction': 0}, 'Batrachomatus': {'total': 5, 'count': 0, 'fraction': 0},
       'Matus': {'total': 4, 'count': 0, 'fraction': 0}}

parser = argparse.ArgumentParser(description="Write file with sample fraction by genus or subfamily")
parser.add_argument("-i", "--input", type=str, help="Text file containing list of taxon names")
parser.add_argument("-t", "--taxonomy", action="store_true", help="Specify level of taxonomy")
args = parser.parse_args()

output = open(f'{args.input}.out', 'w')
output.write('1.0\n')
file = open(args.input)
lines = file.readlines()

# Count taxa in each taxonomic group
for line in lines:
    line = [line.strip()]
    if args.taxonomy == 'subfamily':
        for subfamily, data in sub.items():
            if subfamily in line[0]:
                data['count'] += 1
    if args.taxonomy == 'tribe':
        for tribe, data in tribe.items():
            if tribe in line[0]:
                data['count'] += 1
    if args.taxonomy == 'genus':
        for key, value in subgen.items():
            for v in value:
                if v in line[0]:
                    line[0] = line[0].replace(v, key)
        for genus, data in gen.items():
            if genus in line[0]:
                data['count'] += 1

# Round fraction and write to BAMM file
if args.taxonomy == 'subfamily':
    for subfamily, data in sub.items():
        data['fraction'] = round(data['count'] / data['total'], 3)
if args.taxonomy == 'tribe':
    for tribe, data in tribe.items():
        data['fraction'] = round(data['count'] / data['total'], 3)
if args.taxonomy == 'genus':
    for genus, data in gen.items():
        data['fraction'] = round(data['count'] / data['total'], 3)

# Add fraction to line for BAMM file
for line in lines:
    line = [line.strip()]
    if args.taxonomy == 'subfamily':
        for subfamily, data in sub.items():
            if subfamily in line[0]:
                line = line + [subfamily, str(data['fraction'])]
    if args.taxonomy == 'tribe':
        for tribe, data in tribe.items():
            if tribe in line[0]:
                line = line + [tribe, str(data['fraction'])]
    if args.taxonomy == 'genus':
        for genus, data in gen.items():
            if genus in line[0]:
                line = line + [genus, str(data['fraction'])]
    if len(line) <= 1:
        line = line + ['Outgroup', '1']
    output.write('\t'.join(line) + '\n')
    print('\t'.join(line) + '\n')
