# Input list of species names

sp = {'Agabinae':
          {'total': 426,
           'count': 0,
           'fraction': 0},
      'Colymbetinae':
          {'total': 142,
           'count': 0,
           'fraction': 0},
      'Copelatinae':
          {'total': 791,
           'count': 0,
           'fraction': 0},
      'Coptotominae':
          {'total': 5,
           'count': 0,
           'fraction': 0},
      'Cybistrinae':
          {'total': 129,
           'count': 0,
           'fraction': 0},
      'Dytiscinae':
          {'total': 254,
           'count': 0,
           'fraction': 0},
      'Hydrodytinae':
          {'total': 4,
           'count': 0,
           'fraction': 0},
      'Hydroporinae':
          {'total': 2373,
           'count': 0,
           'fraction': 0},
      'Laccophilinae':
          {'total': 460,
           'count': 0,
           'fraction': 0},
      'Lancetinae':
          {'total': 22,
           'count': 0,
           'fraction': 0},
      'Matinae':
          {'total': 9,
           'count': 0,
           'fraction': 0}}


output = open('test.out', 'w')
output.write('1.0\n')
file = open('test.txt')
lines = file.readlines()
for line in lines:
    line = [line.strip()]
    for subfamily, data in sp.items():
        if subfamily in line[0]:
            data['count'] += 1

for subfamily, data in sp.items():
    data['fraction'] = round(data['count'] / data['total'], 3)

for line in lines:
    line = [line.strip()]
    for subspecies, data in sp.items():
        if subspecies in line[0]:
            line = line + [subspecies, str(data['fraction'])]
    if len(line) <= 1:
        print(line)
        line = line + ['Outgroup', '1']
    output.write('\t'.join(line) + '\n')
