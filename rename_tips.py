import argparse

parser = argparse.ArgumentParser(description="Rename tips of newick tree")
parser.add_argument("-t", "--tree", type=str, help="tree file")
args = parser.parse_args()

tree = (open(args.tree)).read()
tree = tree.split(',')
tips = []
for tip in tree:
    name = tip.split('~', 1)[0]
    data = tip.split("'",2)[2]
    tip = name + "'" + data
    tips.append(tip)
tree = ','.join(tips)
print(tree)
#print(tree)

open(f'rename_{args.tree}', 'w').write(tree)


