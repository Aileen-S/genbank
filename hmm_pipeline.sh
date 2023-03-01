#!/bin/bash
#SBATCH --job-name=hmm_test
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aileen.scott@nhm.ac.uk
#SBATCH -p short
#SBATCH --array=1-10

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
for file in raw/*
do
sed -i 's/ /_/g' $file
done

cat << p

------------
Add Outgorup
------------
p
mv raw/COX1.fasta raw/COX1_without_outgroup.fasta
cat raw/COX1_without_outgroup.fasta outgroup*fasta > raw/COX1.fasta


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
cp -r aa_align aa_align_frametags
for file in aa_align/*
do
   sed -i -E "s/;frame=[0-9]*(;$)?//" $file
done
echo "Frame tags removed from files in aa_align"

cp -r nt_align nt_align_frametags
for file in nt_align/*
do
   sed -i -E "s/;frame=[0-9]*(;$)?//" $file
done
echo "Frame tags removed from files in nt_align"

cat << p

--------
FastTree
--------
p
FastTree -wag aa_align/COX1.fasta > $taxon.tree

cat << p

-------------------------
mPTP Species Delimitation
-------------------------
p
# conda activate mptp?
mkdir mptp
cd mptp

echo "Get minbr value"
outgroup=$(awk 'NR==1' ../outgroup*txt)
# Get minbr value
mptp --minbr_auto alignment --tree_file $taxon.tree --output_file minbr --outgroup $outgroup

# Copy minbr value from minbr.txt
minbr=$(awk 'NR==1' minbr)

echo "Run mPTP"
# Run PTP
mptp --ml --multi --minbr $minbr --tree_file $taxon.tree --output_file mptp --outgroup $outgroup

echo "Run MCMC"
# Run MCMC sampling to determine statistical significance
mptp --mcmc 1000000 --mcmc_sample 1000 --mcmc_log 1000 --mcmc_runs 2 --multi --minbr $minbr --tree_file ../fasttree/$taxon.tree --output_file MCMC
cd ..


