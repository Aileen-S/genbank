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

mkdir $taxon
cd $taxon

cat << p
--------------------
GenBank HMMER search
--------------------
p
#Search GenBank files (downloaded 02/11/23)
mkdir genbank
cd genbank
nhmmer --tblout $taxon.hmmhits ../$taxon.profile.hmm /mnt/shared/projects/nhm/voglerlab/genbank/$taxon.gb

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

python3 ~/scratch/github/genbank/search_gb_file.py -g /mnt/shared/projects/nhm/voglerlab/genbank/$taxon.gb -a $taxon.accessionlist.txt -m -i both
mkdir 1_raw
mv *fasta 1_raw

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
mkdir frame
cd bold

nhmmer --tblout $taxon.hmmhits ../$taxon.profile.hmm /mnt/shared/projects/nhm/voglerlab/bold/*$taxon.fasta

# Save hits as new fasta
python3 ~/scratch/github/genbank/filter_fasta.py -i /home/ascott/voglerlab/bold/*$taxon.fasta -f $taxon.fasta -s $taxon.hmmhits -t hmmer

# Filter fasta to remove dulicate BINs, GenBank sequences and other genes
Rscript ~/scratch/github/voglerlab/bold_cli.R -c /home/ascott/voglerlab/bold/*$taxon.csv -s $taxon.fasta -m $taxon.csv -f $taxon.filter.fasta -g

cat << p

-----------------
Remove Duplicates
-----------------
p
#(maybe do this at AA stage?)

# Merge BOLD and GenBank sequences
cat $taxon.filter.fasta ../genbank/1_raw/COX1.fasta > bold_gb.fasta

# Remove duplicate sequences
echo 'Seaching for duplicates in merged BOLD/GenBank fasta'

python3 ~/scratch/github/genbank/remove_dups.py -i bold_gb.fasta -o nodups.fasta -d duplicates.fasta

cat << p

--------------
Get frame tags
--------------
p
sed -i -E "s/;frame=[0-9]*(;$)?//" nodups.fasta

# Find frame
~/scratch/github/biotools/findframe.py -r ~/voglerlab/6000/$taxon/4_aa_align/COX1.fasta -t 5 < nodups.fasta > frame.fasta
cd ..

cat << p

------------------
Add 6000 sequences
------------------
p
mkdir alignment
cd alignment
mkdir 1_raw
mkdir 2_all

cp ../genbank/1_raw/* 1_raw
cp ../bold/frame.fasta 1_raw/COX1.fasta

# Merge with 6000 fastas for back translation
cd 1_raw
for file in *
do
  cat ~/voglerlab/6000/$taxon/1_raw/${file#*/} $file > ../2_all/${file#*/}
done
cd ..

cat << p

--------------------
Translate to Protein
--------------------
p
mkdir 3_aa
for file in 1_raw/*
do
  echo $file
  ~/scratch/github/biotools/translate.py 5 < $file > 3_aa/${file#*/}
done

# Remove duplicate COX1 amino acid sequences
python3 ~/scratch/github/genbank/remove_dups.py -i 3_aa/COX1.fasta -o nodups.fasta -d duplicates.fasta
cp nodups.fasta 3_aa/COX1.fasta

# Remove duplicates from corresponding raw files
python3 ~/scratch/github/genbank/filter_fasta.py -i 1_raw/COX1.fasta -n 1_raw/COX1.fasta.f -s duplicates.fasta -t fasta
mv 1_raw/COX1.fasta.f 1_raw/COX1.fasta
python3 ~/scratch/github/genbank/filter_fasta.py -i 2_all/COX1.fasta -n 2_all/COX1.fasta.f -s duplicates.fasta -t fasta
mv 2_all/COX1.fasta.f 2_all/COX1.fasta

cat << p

-------------
Align to 6000
-------------
p
mkdir 4_aa_align
for file in 3_aa/*
do
  echo $file
  mafft --add $file --6merpair --maxiterate 1000 --anysymbol --thread 10 ~/voglerlab/6000/$taxon/4_aa_align/${file#*/} >  4_aa_align/${file#*/}
done
echo 'Check alignments, remove any unaligned sequences and proceed to supermatrix pipeline'


