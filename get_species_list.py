# Get list of species (db_xref) from GenBank and remove duplicates
# Each species is a dict
# For each species, each gene is a dict containing all the sequence lenghts
# Find longest, extract and save as fasta

import argparse
import urllib
from Bio import Entrez
from Bio import SeqIO
import textwrap as _textwrap
from collections import defaultdict

# Function definitions

def loadnamevariants():
    conversion = {}
    url = "https://raw.githubusercontent.com/tjcreedy/constants/master/gene_name_variants.txt"
    fullparse = {}
    alltypes = set()
    for line in urllib.request.urlopen(url):
        line = line.decode('utf-8').strip()
        description, variants = line.split(":")
        name, annotype, fullname = description.split(";")
        variants = variants.split(',')
        variants.extend([name, fullname.upper()])

        fullvariants = []
        for v in [name] + variants:
            for g in ['', ' ']:
                v = v.replace(g, '')
                for s in ['', ' GENE', ' '+annotype.upper()]:
                    fullvariants.append(v+s)
                    conversion[v+s] = name

        alltypes.add(annotype)
        fullparse[name] = {'type': annotype, 'variants': fullvariants, 'product': fullname}
    return conversion, alltypes, fullparse


def get_feat_name(feat):
    featname = "unknown"
    nametags = ['gene', 'product', 'label', 'standard_name']  # Search these four keys for gene name
    if any(t in feat.qualifiers.keys() for t in nametags):
        for t in nametags:
            if t in feat.qualifiers.keys():
                featname = feat.qualifiers[t][0].upper()
                break
    return featname


def set_feat_name(feat, name):
    nametags = ['gene', 'product', 'label', 'standard_name']
    if any(t in feat.qualifiers.keys() for t in nametags):
        for t in nametags:
            if t in feat.qualifiers.keys() :
                feat.qualifiers[t][0] = name
    return feat


def search_nuc(term, summaries=False, chunk=10000):
    # Get initial count of responses
    searchhand = Entrez.esearch(db="nucleotide", term=term, retmax=0)
    searchrec = Entrez.read(searchhand)
    count = int(searchrec["Count"])

    # Yield
    for start in range(0, count, chunk):
        # Search and get GB IDs
        searchhand = Entrez.esearch(db="nucleotide", term=term, retstart=start, retmax=chunk)
        searchrec = Entrez.read(searchhand)
        gbids = searchrec['IdList']
        # Yield only GBIDs if no summaries desired
        if not summaries:
            yield gbids
        else:
            # Retrieve summaries and yield both otherwise
            sumhand = Entrez.esummary(db="nucleotide", id=','.join(gbids))
            sumrec = Entrez.read(sumhand)
            yield gbids, sumrec


# This reclasses the argparse.HelpFormatter object to have newlines in the help text for paragraphs
class MultilineFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(' ', text).strip()
        paragraphs = text.split('|n ')
        multiline_text = ''
        for paragraph in paragraphs:
            formatted_paragraph = _textwrap.fill(paragraph, width, initial_indent=indent,
                                                 subsequent_indent=indent
                                                 ) + '\n\n'
            multiline_text = multiline_text + formatted_paragraph
        return multiline_text

# Argument parser
parser = argparse.ArgumentParser(description="Search GenBank, retrive gene sequences and save as fasta.", formatter_class=MultilineFormatter)

parser.add_argument("-t", "--taxon", type=str, help="Taxon of interest: must be specified")
parser.add_argument("-g", "--gene", type=str, help="Gene(s) of interest: if not specified, all genes will be retrieved")
parser.add_argument("-e", "--email", type=str, help="Your email registered with NCBI")
parser.add_argument("-n", "--nuclear", action="store_true", help="Search for nuclear as well as mitochondrial genes.")
# parser.add_argument("-l", "--length", type=str)


# Start the actual script

args = parser.parse_args()         # Process input args from command line
# args = parser.parse_args('-t Dytiscidae -e '.split(' ')) # This is how I step through the script interactively

# Get name variants
nameconvert, types, namevariants = loadnamevariants()
# So issue 1 was that you had an old version of the loadnamevariants function, sorry about that.
# I've updated this here - you can just ignore the types and namevariants objects, you (probably)
# won't need them

# Set up for unrecognised genes
unrecgenes = defaultdict(list)


Entrez.email = args.email

# Generate search term to get all sequences in the search taxonomy
# - if -n option not used, then include "mitochondrial" in search term.
basesearch = f"(\"{args.taxon}\"[Organism] OR \"{args.taxon}\"[All Fields])"\ 
             f"{'' if args.nuclear else ' AND mitochondrion[filter]'}"

# Retrieve all taxids that represent tips of the NCBI Taxonomy tree
# Make the search generator
searchgen = search_nuc(term=basesearch, summaries=True, chunk=5000)
taxids = set()
#i = 0
for gbids, summaries in searchgen:
    # gbids, summaries = next(searchgen)
    #i += 1
    taxa = set(int(s['TaxId']) for s in summaries)
    taxids.update(taxa)
    #print(f"iteration={i}, returns={len(gbids)}, first gbid={gbids[0]}, first summary accession={summaries[0]['Caption']}, taxids in this iteration={len(taxa)}, total taxids={len(taxids)}")

# Some of these will be subspecies.
# You need to search them in NCBI Taxonomy to weed out the subspecies and generate a list of latin biomials.
# Then iterate through each of these binomials (not taxids as initially thought) to download the sequences etc


print(str(len(taxids)) + " unique species saved")                       # Print total taxon ids

species = {}
for tax in taxids:
    handle = Entrez.esearch(db="nucleotide", term=tax)                # Search for all records for each taxon id
    # Issue 2 is that you're not narrowing down your search here to include only mitochondrial
    # sequences. Hence you'll get some weird names in your unrecgenes like ARK, H3 or H4 - these
    # are nuclear genes. Now of course at some point you might want these, but the namevariants
    # file is not set up for these so they will be ignored for now
    record = Entrez.read(handle)
    accs   = record["IdList"]                                         # Get accessions
    accstr = ",".join(accs)                                           # Join into string for efetch
    handle = Entrez.efetch(db="nucleotide", id=accstr, rettype="gb", retmode="text")  # Get GenBanks
    record = SeqIO.parse(handle, "gb")
    for rec in record:
        for feature in rec.features:
            type = feature.type                                 # Retrieve the feature type, might be useful later
            # Issue 3 is that this is useful now! If you inspect any genbank file, you'll see that
            # most genes have both a gene feature and a CDS/tRNA/rRNA feature. This makes sense for
            # nuclear sequences but is meaningless duplication in mitogenomes. There'll also be
            # various other features, the standard one being a 'source' annotation which contains
            # more metadata (I don't really know why either...). So, to remove duplicates and
            # ignore other types of annotation, I would add a filter here to only process features
            # that are CDS or rRNA, and just skip any others. Alternatively, you could just look at
            # gene annotations, although you'd need to filter out the tRNAs at some point...
            if type not in ('CDS', 'rRNA'):
                continue
            # continue = skip the rest of the current iteration of this loop, i.e. in this case, go
            # on to the next feature in rec.features
            name = get_feat_name(feature)                       # Use function to search for gene names
            if name in nameconvert:
                stdname = nameconvert[name]                    # If gene name in namevariants, convert to standard name
                # You might want to add a catch here to only include features if they are in a list
                # of required genes that you supply.
                sequence = rec[feature.location.start:feature.location.end]
                output = [rec.name, type, len(sequence)]        # Save relevent info in list (Add sequence when working)
                if tax in species:                              # If taxon ID in dict
                    if stdname in species[tax]:                 # If gene in dict for that taxon ID
                        species[tax][stdname].append(output)    # Add gene info list to dict
                    else:
                        species[tax] = {stdname: [output]}
                else:
                    species[tax] = {stdname: [output]}      # Otherwise add to dict with new key

            else:
                unrecgenes[rec.name].append(name)            # If gene name not in namevarants, save to list to check later
                # I changed the key here to rec.name so it makes a little more sense to me

print("\nSpecies Dict")
print(species)
print("\nUnrecognised Genes")
print(unrecgenes)


# for k,v in Species.items():

# Add mitochondrion to term
# Max sequence length 20000[SLEN]