#!/bin/bash
#SBATCH --job-name=hmm_test
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aileen.scott@nhm.ac.uk
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-4


# First get database to search, and profile to search with.
# Save profile as taxon.profile.fasta

# In directory named for each taxon as in config file:
#   Save outgroup COX1 sequences as "outgroup.fasta"
#   Save outgroup taxa names seperated by a comma

# If using slurm, get taxon names from config file
taxon=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' taxa.config)

cd $taxon

cat << p
--------------------
GenBank HMMER search
--------------------
p
#Search GenBank files (downloaded 02/11/23)
mkdir genbank
cd genbank
nhmmer --tblout $taxon.hmmhits ../../profile.hmm /mnt/shared/projects/nhm/voglerlab/genbank/$taxon.gb

cat << p
------------------------------
Get list of GenBank accessions
------------------------------
p
python3 ~/scratch/github/genbank/hmm_get_accessions.py -i $taxon.hmmhits -o $taxon.accessionlist.txt

cat << p

--------------------------------------------
Get sequences and metadata from GenBank file
--------------------------------------------
p

python3 ~/scratch/github/genbank/search_gb_file.py -g /mnt/shared/projects/nhm/voglerlab/genbank/$taxon.gb -a ../$taxon.accessionlist.txt -m -i both
mv metadata.csv ..

# Remove spaces
for file in 1_raw/*
do
sed -i 's/ /_/g' $file
done

cd ..

cat << p

-----------------
BOLD HMMER search
-----------------
p

mkdir bold
cd bold

nhmmer --tblout $taxon.hmmhits ../../profile.hmm /mnt/shared/projects/nhm/voglerlab/bold/*$taxon.fasta

# Save hits as new fasta
python3 ~/scratch/github/genbank/filter_fasta.py -i /home/ascott/voglerlab/bold/*$taxon.fasta -o $taxon.fasta -hm $taxon.hmmhits

# Filter fasta to remove dulicate BINs, GenBank sequences and other genes
Rscript ~/scratch/github/voglerlab/bold_cli.R -c /home/ascott/voglerlab/bold/*$taxon.csv -s $taxon.fasta -m $taxon.csv -f $taxon.filter.fasta -g

cat << p

---------------------------------
Get frame tags for BOLD sequences
---------------------------------
p

mafft --add $taxon.filter.fasta --maxiterate 1000 --thread 10 /mnt/shared/projects/nhm/voglerlab/6000/$taxon/5_nt_align/COX1.fasta > $taxon.align.fasta
~/scratch/github/biotools/findframe.py -a $taxon.align.fasta < $taxon.filter.fasta > COX1.fasta
cd ..

cat << p

--------------
Join sequences
--------------
p

mkdir alignment
cd alignment

# Merge with 6000 fastas for back translation
mkdir 1_raw
for file in /mnt/shared/projects/nhm/voglerlab/6000/$taxon/1_raw/*
do
  cat $file ../genbank/${file#*/} ../bold/${file#*/} > 1_raw/${file#*/}
done

# Merge bold and genbank for pipeline
mkdir 2_merge
for file in 1_raw/*
do
  cat ../genbank/${file#*/} ../bold/${file#*/} > 2_merge/${file#*/}
done

cat << p

--------------------
Translate to Protein
--------------------
p
mkdir 3_aa
for file in 2_merge/*
do
  echo $file
  ~/scratch/github/biotools/translate.py 5 < $file > 3_aa/${file#*/}
done

cat << p

-------------
Align to 6000
-------------
p
mkdir 4_aa_align
for file in 3_aa/*
do
  echo $file
  mafft --add $file --6merpair --maxiterate 1000 --anysymbol --thread 10 /mnt/shared/projects/nhm/voglerlab/6000/$taxon/4_aa_align/${file#*/} >  4_aa_align/${file#*/}
done

cat << p

---------------------------------
Translate to Nucleotide Alignment
---------------------------------
p
mkdir 5_nt_align
for file in 4_aa_align/*
do
  echo $file
  ~/scratch/github/biotools/backtranslate.py -i $file 2_merge/${file#*/} 5 > 5_nt_align/${file#*/}
done

echo 'Check alignments and proceed to PTP pipeline'