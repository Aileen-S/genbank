#!/bin/bash
#SBATCH --job-name=ptp
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

----------------
Write constraint
----------------
p
cd supermatrix

# Constraint specifies outgroup for mPTP
rm outgroup*
cp /home/ascott/voglerlab/6000/$taxon/outgroup6000.fasta  .
sed -i -E "s/;frame=[0-9]*(;$)?//" outgroup6000.fasta
python3 ~/scratch/github/genbank/filter_fasta.py -i 1_aa_supermatrix.fasta -f outgroup.fasta -n ingroup.fasta -s outgroup6000.fasta -t fasta

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
FastTree -gtr -nt -constraints ../supermatrix/fasttree_constraint.txt < ../supermatrix/2_nt_supermatrix.fasta > $taxon.nt.tree
FastTree -wag -constraints ../supermatrix/fasttree_constraint.txt < ../supermatrix/1_aa_supermatrix.fasta > $taxon.aa.tree

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


cat << p

-------
IQ_TREE
-------
p

mkdir iqtree
cd iqtree

# Run condon tree with MFP+MERGE
iqtree -s ../supermatrix/2_nt_supermatrix.fasta --prefix codon -m MFP+MERGE -sp ../supermatrix/partitions_codon123.txt -g ~/voglerlab/6000/$taxon/new.tree
iqtree -s ../supermatrix/1_aa_supermatrix.fasta --prefix aa -m MFP+MERGE -sp ../supermatrix/partitions_aa.txt -g ~/voglerlab/6000/$taxon/new.tree

cd ..