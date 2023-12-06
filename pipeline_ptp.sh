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
cd alignment

cat << p

---------------
Filter raw COX1
---------------
p
# Assumes 4_aa_align has been checked and unaligned sequences have been removed

# Filter 2_all COX1 file for backtranslation
python3 ~/scratch/github/genbank/filter_fasta.py -i 2_all/COX1.fasta -f 2_all/COX1.f -s 4_aa_align/COX1.fasta -t fasta
mv 2_all/COX1.f 2_all/COX1.fasta

cat << p

---------------------------------
Translate to Nucleotide Alignment
---------------------------------
p
mkdir 5_nt_align

for file in 4_aa_align/*
do
  echo $file
  ~/scratch/github/biotools/backtranslate.py -s -i $file 2_all/${file#*/} 5 > 5_nt_align/${file#*/}
done

cat << p

-----------------
Remove Frame Tags
-----------------
p

for file in 4_aa_align/*
do
  echo $file
  sed -i -E "s/;frame=[0-9]*(;$)?//" $file
done
echo "Frame tags removed"

cat << p

-----------------------------------------
Get alignment files ready for supermatrix
-----------------------------------------
p

mkdir ../supermatrix
mkdir ../supermatrix/aa
mkdir ../supermatrix/nt

for file in 4_aa_align/*
do
  python3 ~/scratch/github/genbank/fasta_remove_gbids.py -i $file -o ../supermatrix/aa/${file#*/}
done

for file in 5_nt_align/*
do
  python3 ~/scratch/github/genbank/fasta_remove_gbids.py -i $file -o ../supermatrix/nt/${file#*/}
done

# Filter out any other genes for taxa removed from COX1
for file in 2_all/*
do
  echo $file
  python3 ~/scratch/github/genbank/filter_fasta.py --input $file --found 2_all/${file#*/}.f --search 4_aa_align/COX1.fasta --type fasta
  mv 2_all/${file#*/}.f 2_all/${file#*/}
done

cd ../supermatrix

for file in aa/*
do
  echo $file
  python3 ~/scratch/github/genbank/filter_fasta.py -i $file -f aa/${file#*/}.f -s ../alignment/4_aa_align/COX1.fasta -t fasta
  mv aa/${file#*/}.f aa/${file#*/}
done

for file in nt/*
do
  echo $file
  python3 ~/scratch/github/genbank/filter_fasta.py -i $file -f nt/${file#*/}.f -s ../alignment/4_aa_align/COX1.fasta -t fasta
  mv nt/${file#*/}.f nt/${file#*/}
done

cat << p

---------------------------------------
Build Supermatrices and Partition Files
---------------------------------------
p

~/scratch/github/catfasta2phyml/catfasta2phyml.pl -c -fasta aa/* > 1_aa_supermatrix.fasta 2> 1_aa_partitions.txt
~/scratch/github/catfasta2phyml/catfasta2phyml.pl -c -fasta nt/* > 2_nt_supermatrix.fasta 2> 2_nt_partitions.txt

python3 ~/scratch/github/genbank/partitions.py -i 1_aa_partitions.txt -t aa -o raxml
python3 ~/scratch/github/genbank/partitions.py -i 2_nt_partitions.txt -t nt -o raxml


# Change partitions to match a previous run (harder with codon model. Postpone!

#sed -i -E "s/aa\///g" 1_aa_partitions.txt
#sed -i -E "s/nt\///g" 2_nt_partitions.txt
#sed -i -E "s/.fasta//g" *_partitions.txt
#sed -i -E "s/.fasta//g" ../iqtree/*.best_scheme.nex

#~/scratch/github/phylostuff/partitioner.py -p ../iqtree/nttree_MF.best_scheme.nex < 2_nt_partitions.txt > partitions_gene.nex
#~/scratch/github/phylostuff/partitioner.py -p ../iqtree/aa_MF.best_scheme.nex < 1_aa_partitions.txt > partitions_aa.nex
cat << p

----------------
Write constraint
----------------
p

# Constraint specifies outgroup for mPTP
rm outgroup*
cp /home/ascott/voglerlab/6000/$taxon/outgroup*  .
sed -i -E "s/;frame=[0-9]*(;$)?//" outgroup*
python3 ~/scratch/github/genbank/filter_fasta.py -i 1_aa_supermatrix.fasta -f outgroup.fasta -n ingroup.fasta -s outgroup* -t fasta

# Outputs raxml_constraint.txt, fasttree_constraint.txt and outgroup.txt
python3 ~/scratch/github/genbank/write_constraint.py -i ingroup.fasta -o outgroup.fasta

cd ..

cat << p

--------
FastTree
--------
p

mkdir mptp
cd mptp

# Run FastTree
FastTree -gtr -nt -constraints ../supermatrix/fasttree_constraint.txt < ../supermatrix/2_nt_supermatrix.fasta > $taxon.tree
#FastTree -wag -constraints ../supermatrix/fasttree_constraint.txt < ../supermatrix/1_aa_supermatrix.fasta > $taxon.tree



cat << p

-------
IQ_TREE
-------
p

mkdir iqtree
cd iqtree

# Run condon tree with MFP+MERGE
iqtree -s ../supermatrix/2_nt_supermatrix.fasta --prefix codon -m MFP+MERGE -sp ../supermatrix/partitions_codon123.txt
iqtree -s ../supermatrix/1_aa_supermatrix.fasta --prefix codon -m MFP+MERGE -sp ../supermatrix/partitions_aa.txt

cd ..

cat << p

-------------------------
mPTP Species Delimitation
-------------------------
p

# Check tree and make outgroup file

# Get outgroup
outgroup=$(awk 'NR==1' ../supermatrix/outgroup.txt)

# Calculate minbr value
echo "Getting minbr value"
mptp --minbr_auto ../supermatrix/2_nt_supermatrix.fasta --tree_file $taxon.tree --output_file minbr --outgroup $outgroup

# Copy minbr value from slurm output
minbr=$(awk 'END {print $NF}' ../../slurm*$SLURM_ARRAY_TASK_ID.out)
echo "minbar value = " $minbr

# Run PTP
echo "
Running single mPTP"
mptp --ml --single --minbr  $minbr --tree_file $taxon.tree --output_file mptp_single --outgroup $outgroup

# Run MCMC sampling to determine statistical significance
#echo "
#Running MCMC"
#mptp --mcmc 1000000 --mcmc_sample 1000 --mcmc_log 1000 --mcmc_runs 2 --multi --minbr $minbr --tree_file $taxon.tree --output_file mcmc
#cd ..


