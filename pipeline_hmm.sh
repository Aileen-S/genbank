#!/bin/bash
#SBATCH --job-name=hmm_test
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aileen.scott@nhm.ac.uk
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4
#SBATCH --array=10


# First get database to search, and profile to search with.
# Save profile as taxon.profile.fasta

# In directory named for each taxon as in config file:
# Save outgroup COX1 sequences as "outgroup.fasta"
# Save outgroup taxa names seperated by a comma

# If using slurm, get taxon names from config file
taxon=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' adephaga.txt)

cd $taxon

cat << p

------------
HMMER search
------------
p
# Add -A $taxonout.sth for sequence file output
nhmmer --tformat fasta --tblout $taxon.hmmhits.txt ../../beetle_profile.hmm /home/ascott/voglerlab/barcodes/fastas/Adephaga/$taxon.fasta

cat << p

------------------------------
Get list of GenBank Accessions
------------------------------
p
python3 ~/scratch/github/genbank/hmm_get_accessions.py -i $taxon.hmmhits.txt -o $taxon.accessionlist.txt

cat << p

-----------------------------
Search GenBank for Accessions
-----------------------------
p
mkdir raw
cd raw
python3 ~/scratch/github/genbank/get_concat_recs.py -f ../$taxon.accessionlist.txt -r gbid -i both -e mixedupvoyage@gmail.com -m
mv metadata.csv ..
cd ..

# Remove spaces
for file in raw/*
do
sed -i 's/ /_/g' $file
done


cat << p

----------------
Write Constraint
----------------
p
mkdir trees
cd trees

python3 ~/scratch/github/genbank/write_constraint.py -i ../raw/COX1.fasta -o ../outgroup*fasta
#perl ~/scratch/github/fasttree_constraint.pl < raxml_constraint.txt > constraint_alignment.txt

# Remove Frame Tage
#sed -i -E "s/;frame=[0-9]*(;$)?//" raxml_constraint.txt
sed -i -E "s/;frame=[0-9]*(;$)?//" fasttree_constraint.txt
sed -i "s/;frame=[0-9]//g" raxml_constraint.txt

cd ..

cat << p

------------
Add Outgroup
------------
p
mv raw/COX1.fasta COX1raw.fasta
cat COX1raw.fasta outgroup*fasta > raw/COX1.fasta


cat << p

--------------------
Translate to Protein
--------------------
p
mkdir aa
for file in raw/*
do
  echo $file
  ~/scratch/github/biotools/translate.py 5 < $file > aa/${file#*/}
done

cat << p

----------------
Align to Profile
----------------
p
mkdir aa_align
for file in aa/*
do
  echo $file
  mafft --add $file --6merpair --maxiterate 1000 --anysymbol --thread 10 ~/scratch/profiles/0_AA_profiles/${file#*/} >  aa_align/${file#*/}
done


cat << p

--------------
Remove Profile
--------------
p
cd aa_align
for file in *
do
  echo $file
  mv $file ${file}_bak
  perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' ${file}_bak | grep -A 1 -P '^>(?!PROFILE::).*' > $file
  rm ${file}_bak
done
cd ..

cat << p

---------------------------------
Translate to Nucleotide Alignment
---------------------------------
p
mkdir nt_align
for file in aa_align/*
do
  echo $file
  ~/scratch/github/biotools/backtranslate.py -i $file raw/${file#*/} 5 > nt_align/${file#*/}
done

cat << p

-----------------
Remove Frame Tags
-----------------
p
#cp -r aa_align aa_align_frametags
for file in aa_align/*
do
   sed -i -E "s/;frame=[0-9]*(;$)?//" $file
done
echo "Frame tags removed from files in aa_align"

#cp -r nt_align nt_align_frametags
for file in nt_align/*
do
   sed -i -E "s/;frame=[0-9]*(;$)?//" $file
done
echo "Frame tags removed from files in nt_align"

cat << p

---------------------------------------
Build Supermatrices and Partition Files
---------------------------------------
p

mkdir trees
cd trees

~/scratch/github/catfasta2phyml/catfasta2phyml.pl -c -fasta aa_align/* > 1_aa_supermatrix.fasta 2> 1_aa_partitions.txt
~/scratch/github/catfasta2phyml/catfasta2phyml.pl -c -fasta nt_align/* > 2_nt_supermatrix.fasta 2> 2_nt_partitions.txt

cat << p

-----
RAxML
-----
p

# Get random number seed for RAxML and save to slurm output
seed=$RANDOM
echo 'RAxML random number seed = '$seed

# Run RaxML
raxml-ng --msa ../2_nt_supermatrix.fasta --model GTR+G --prefix $taxon --threads 2 --seed $seed --tree-constraint raxml_constraint.txt


cat << p

--------
FastTree
--------
p
# Run FastTree
FastTree -gtr -nt -constraints fasttree_constraint.txt < ../2_nt_supermatrix.fasta > $taxon.tree
cd ..

cat << p

------------------------
PTP Species Delimitation
------------------------
p
mkdir ptp
cd ptp

seed=$RANDOM
echo 'PTP random number seed = '$seed

bPTP.py -t ../trees/$taxon.raxml.bestTree -o bPTP -s $seed

cd ..


cat << p

-------------------------
mPTP Species Delimitation
-------------------------
p
mkdir mptp
cd mptp

# Get outgroup
outgroup=$(awk 'NR==1' ../outgroup*txt)

# Calculate minbr value
echo "Getting minbr value"
mptp --minbr_auto ../nt_align/COX1.fasta --tree_file ../trees/$taxon.raxml.bestTree --output_file minbr --outgroup $outgroup

# Copy minbr value from slurm output
minbr=$(awk 'END {print $NF}' ../../slurm*$SLURM_ARRAY_TASK_ID.out)
echo "minbar value = " $minbr

# Run PTP
echo "
Running mPTP"
mptp --ml --multi --minbr  $minbr --tree_file ../trees/$taxon.raxml.bestTree --output_file mptp --outgroup $outgroup

# Run MCMC sampling to determine statistical significance
echo "
Running MCMC"
mptp --mcmc 1000000 --mcmc_sample 1000 --mcmc_log 1000 --mcmc_runs 2 --multi --minbr $minbr --tree_file ../trees/$taxon.raxml.bestTree --output_file mcmc
cd ..


