#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Description of what the script does"""

# Imports

import argparse
import textwrap as _textwrap

from Bio import SeqIO

# Class definitions

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

# Global variables

# Function definitions

def loadnamevariants():
    output = {}
    url = "https://raw.githubusercontent.com/tjcreedy/genenames/master/gene_name_variants.txt"
    for line in urllib.request.urlopen(url):
        line = line.decode('utf-8').strip()
        name = line.split(";")[0]
        annotype = line.split(":")[0].split(";")[1]
        variants = line.split(":")[1].split(",")
        for v in variants:
            for g in ['', ' ']:
                v = v.replace(g, '')
                for s in ['',' GENE', ' '+annotype.upper()]:
                    output[v+s] = name
    return(output)

def getcliargs(arglist=None):
    # Initialise the parser object and give a description
    parser = argparse.ArgumentParser(description="""
        Here describe what the code does. I often find it useful to draft this before I even write
        the code, because it helps me think things through from the end-user perspective (but it 
        may be that a little experience is needed before you can think like this and that's OK!)
        |n
        Separate paragraphs with the symbol above. It doesn't have to be on its own line but I find
        it nicer
        """, formatter_class=MultilineFormatter)

    # Add individual argument specifications
    parser.add_argument("-f", "--flag", action='store_true',
                        help="this argument stores True if used and False otherwise")
    parser.add_argument("-p", "--filepath", type=str, metavar="PATH", required=True,
                        help="this argument records the path to a file, it's required")
    parser.add_argument("-n", "--number", type=int, metavar="N", default=2,
                        help="this argument records an integer, it must be 1, 2, 3 or 4. If not "
                             "supplied, 2 will be used",
                        choices=[1,2,3,4])

    # Parse the arguments from the function call if specified, otherwise from the command line
    args = parser.parse_args(arglist) if arglist else parser.parse_args()

    # Do some checking of the inputs
    if args.filepath is not "what you want"
        parser.error(f"{args.filpath} is not what I want!")

    # If the arguments are all OK, output them
    return args

# Start the actual script
if __name__ == "__main__":

    # Get the arguments
    args = getcliargs() # Try to read from the command line, won't work interactively!
    args = getcliargs('-f -p mypath.txt -n 4'.split(' ')) # Read from a string, good for testing




    # First step: search genbank
        # Retrieve the gene name variants
        varianttoname = loadnamevariants()
        nametovariant = # Do something to convert to dict of 13 lists variants for each name

    #say taxsearch = some list of taxa to search genbank for

    allaccessions = []
    for taxon in taxsearch:
        for gene in nametovariant.keys():

            # Run your esearch - btw check out Bio.Entrez.esearch
            searchresults = esearch(...)
            allaccessions.append(searchresults)

    allaccessions = set(allaccessions) # find the unique by converting to a different sort of array

    # Retrieve the genbank data

    # here you could steal my get_genbanks code
    # or you could use Bio.Entrez.efetch - checkout my get_NCBI_taxonomy script for this, you may
    # need to have some more authentication info but see the functions in that rscripts

    efetch(allaccessions, genbankpath)

    # Parse the genbank file
    # The key thing here is that we're only looking at each entry in the genbank file once, to
    # save memory and time


    # Set up a dict for the species data
    speciesdata = {}

    for record in SeqIO.parse(genbankpath, "genbank"):
        # Note - to get a record for code development and testing run record = next(SeqIO.parse(genbankpath, "genbank"))
        # Get the key metadata, i.e. the accession number and species name

        # Check to see if you already have an entry for this species, if not create an entry
        # comprising an empty dict

        # Iterate through the sequence features (aka sequence annotations)
        for feat in record.features:
            # To get a feature for code development and testing run feat = record.features[0]
            # Skip this feature if it's not a CDS (hint: look up the 'continue' command

            # Convert the feature id into the standard names
            # Not, you might want to use try except here to throw an error to the user if the
            # feature has a name that's not in my variants table

            # Get the feature length

            # Create a 3-entry list of the accession number, feature id and length
            # Check to see if you already have an entry for this gene for this species (hint: speciesdata[species][gene]):
                # if it's absent, create an entry comprising a list of one item, that item being the 3-entry list you already made
                # if it's present, append your 3-entry list

    # Filter the data for the ?longest? sequence for each gene in each species

    # Set up a dict of genes from accessions to use
    selecteddata = {}

    for species, data in speciesdata.items():
        for gene, seqs in data.items():

            # Get all the lengths - hint, use list comprehension or a loop to find the nth value in each list in a list of lists

            # Get the maximum length (or whichever length you want, perhaps the method you use is
            # determined by a command line argument?

            # Get the entry or entries from seqs with the selected length - hint, use list comprehension or a loop

            # What will you do if multiple entries have the same target length?

            # Check to see if this accession has an entry in selecteddata already
                # If it's absent, create a dict comprising a dict of genes, e.g.
                selecteddata[accession] = {variant: standardname}
                # If it's present, you just need to add the gene to the dict of genes, e.g.
                selectedata[accession][variant] = standardname

    # Finally, write out your genes
        # Set up a dict of file handles to write to (none of this with open shit here)
    fh = {}
    for gene in nametovariant.keys():
        fh[gene] = open(path, 'w')

        # Go through your genbank file again
    for record in SeqIO.parse(genbankpath, "genbank"):

        # Skip the accession if not in your selected data

        # Get the species name

        # Go through the features
        for feat in record.features:
            # Skip the feature if the id is not in selecteddata[accession][genes].keys()

            # Extract the sequence


            # Look up what standard gene name it is

            # Write the sequence out to the appropriate handle
            fh[gene].write(f">{selecteddata[accession][species]}\nsequence\n")

    # Close your files
    for gene, handle in fh.items():
        handle.close()

    # You're done!
