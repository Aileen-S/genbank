#!/bin/bash
#SBATCH --job-name=blast_test
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aileen.scott@nhm.ac.uk
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-4

# Setup:
# Save taxa in taxa.config config file
# Save search profile as profile.fasta
# In directory named for each taxon as in config file:
#   Save outgroup COX1 sequences as "outgroup.fasta"
#   Save outgroup taxa names seperated by a comma

# Get taxon names from config file
taxon=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' taxa.config)
txid=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' taxa.config)

cd $taxon


cat << p

---------------------
Search BLAST Database
---------------------
p

# Direct to BLAST database
export BLASTDB=/mnt/shared/apps/databases/ncbi/

# Get TXID list
get_species_taxids.sh -t $txid > $taxon.txids

# Search BLAST nucleotide database using profile, and limit search by txid. Output list of accessions.
blastn -db nt -query ../profile.fasta -taxidlist $taxon.txids -out $taxon.blast -max_target_seqs 100000 -outfmt '6 sacc'


cat << p

-----------------------------
Search GenBank for Accessions
-----------------------------
p
mkdir 1_raw
cd 1_raw
python3 ~/scratch/github/genbank/get_concat_recs.py -f ../$taxon.blast -r gbid -i both -e aileen.scott@nhm.ac.uk -m
mv metadata.csv ..
cd ..

# Remove spaces
for file in 1_raw/*
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

python3 ~/scratch/github/genbank/write_constraint.py -i ../1_raw/COX1.fasta -o ../outgroup/COX1.fasta
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
mkdir 2_merge
for file in 1_raw/*
do
  cat $file outgroup/${file#*/} > 2_merge/${file#*/}
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

----------------
Align to Profile
----------------
p
mkdir 4_aa_align
for file in 3_aa/*
do
  echo $file
  mafft --add $file --6merpair --maxiterate 1000 --anysymbol --thread 10 ~/scratch/profiles/0_AA_profiles/${file#*/} >  4_aa_align/${file#*/}
done


cat << p

--------------
Remove Profile
--------------
p
cd 4_aa_align
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
mkdir 5_nt_align
for file in 4_aa_align/*
do
  echo $file
  ~/scratch/github/biotools/backtranslate.py -i $file 2_merge/${file#*/} 5 > 5_nt_align/${file#*/}
done

cat << p

-----------------
Remove Frame Tags
-----------------
p
#cp -r 4_aa_align aa_align_frametags
for file in 4_aa_align/*
do
   sed -i -E "s/;frame=[0-9]*(;$)?//" $file
done
echo "Frame tags removed from files in 4_aa_align"

#cp -r 5_nt_align nt_align_frametags
for file in 5_nt_align/*
do
   sed -i -E "s/;frame=[0-9]*(;$)?//" $file
done
echo "Frame tags removed from files in 5_nt_align"

cat << p

---------------------------------------
Build Supermatrices and Partition Files
---------------------------------------
p

mkdir supermatrix
mkdir supermatrix/aa
mkdir supermatrix/nt

for file in 4_aa_align/*
do
  python3 ~/scratch/github/genbank/fasta_remove_gbids.py -i $file -o supermatrix/aa/${file#*/}
done

for file in 5_nt_align/*
do
  python3 ~/scratch/github/genbank/fasta_remove_gbids.py -i $file -o supermatrix/nt/${file#*/}
done

cd supermatrix
~/scratch/github/catfasta2phyml/catfasta2phyml.pl -c -fasta aa/* > 1_aa_supermatrix.fasta 2> 1_aa_partitions.txt
~/scratch/github/catfasta2phyml/catfasta2phyml.pl -c -fasta nt/* > 2_nt_supermatrix.fasta 2> 2_nt_partitions.txt

python3 ~/scratch/github/genbank/partitions.py -i 1_aa_partitions.txt -t aa -o raxml
python3 ~/scratch/github/genbank/partitions.py -i 2_nt_partitions.txt -t nt -o raxml
cd ..


# Change partitions to match a previous run (harder with codon model. Postpone!

sed -i -E "s/aa\///g" 1_aa_partitions.txt
sed -i -E "s/nt\///g" 2_nt_partitions.txt
sed -i -E "s/.fasta//g" *_partitions.txt
sed -i -E "s/.fasta//g" ../iqtree/*.best_scheme.nex

~/scratch/github/phylostuff/partitioner.py -p ../iqtree/nttree_MF.best_scheme.nex < 2_nt_partitions.txt > partitions_gene.nex
~/scratch/github/phylostuff/partitioner.py -p ../iqtree/aa_MF.best_scheme.nex < 1_aa_partitions.txt > partitions_aa.nex


cat << p

-----
RAxML
-----
p

cd supermatrix
python3 ~/scratch/github/genbank/partitions_raxml.py -i partitions_gene.nex
cd ..

cd trees

# Get random number seed for RAxML and save to slurm output
seed=$RANDOM
echo 'RAxML random number seed = '$seed

# Run RaxML
raxml-ng --msa ../supermatrix/2_nt_multi.fasta --model ../supermatrix/partitions_gene.nex.raxml --prefix $taxon --threads 2 --seed $seed #--tree-constraint raxml_constraint.txt


cat << p

--------
FastTree
--------
p
# Run FastTree
FastTree -gtr -nt -constraints fasttree_constraint.txt < ../supermatrix/2_nt_supermatrix.fasta > $taxon.tree
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

bPTP.py -t ../trees/$taxon.tree -o bPTP -s $seed

cd ..


cat << p

-------------------------
mPTP Species Delimitation
-------------------------
p
mkdir mptp
cd mptp

# Get outgroup
outgroup=$(awk 'NR==1' ../outgroup.txt)

# Calculate minbr value
echo "Getting minbr value"
mptp --minbr_auto ../5_nt_align/COX1.fasta --tree_file ../trees/$taxon.raxml.bestTree --output_file minbr --outgroup $outgroup

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


