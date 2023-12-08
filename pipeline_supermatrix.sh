#!/bin/bash
#SBATCH --job-name=hmm_test
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aileen.scott@nhm.ac.uk
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-4

# Follows on from pipeline_hmm.sh
# Assumes 4_aa_align has been checked and unaligned sequences have been removed from 4_aa_align/COX1.fasta

taxon=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' taxa.config)

cd $taxon
cd alignment

cat << p

----------------------
Filter sequences files
----------------------
p

# Copy edited COX1 file as filter and remove frame tags
cp 4_aa_align/COX1.fasta filter.fasta
sed -i -E "s/;frame=[0-9]*(;$)?//" filter.fasta

# Remove any taxa missing from COX1
for file in 1_raw/*
do
  echo $file
  python3 ~/scratch/github/genbank/filter_fasta.py -i $file -f 1_raw/${file#*/}.f -s filter.fasta -t fasta
  #mv 1_raw/${file#*/}.f 1_raw/${file#*/}
done

for file in 2_all/*
do
  echo $file
  python3 ~/scratch/github/genbank/filter_fasta.py -i $file -f 2_all/${file#*/}.f -s filter.fasta -t fasta
  mv 2_all/${file#*/}.f 2_all/${file#*/}
done

for file in 3_aa/*
do
  echo $file
  python3 ~/scratch/github/genbank/filter_fasta.py -i $file -f 3_aa/${file#*/}.f -s filter.fasta -t fasta
  mv 3_aa/${file#*/}.f 3_aa/${file#*/}
done

for file in 4_aa_align/*
do
  echo $file
  python3 ~/scratch/github/genbank/filter_fasta.py -i $file -f 4_aa_align/${file#*/}.f -s filter.fasta -t fasta
  mv 4_aa_align/${file#*/}.f 4_aa_align/${file#*/}
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
cd ..

cat << p

---------------------------------------
Build Supermatrices and Partition Files
---------------------------------------
p

cd supermatrix
~/scratch/github/catfasta2phyml/catfasta2phyml.pl -c -fasta aa/* > 1_aa_supermatrix.fasta 2> 1_aa_partitions.txt
~/scratch/github/catfasta2phyml/catfasta2phyml.pl -c -fasta nt/* > 2_nt_supermatrix.fasta 2> 2_nt_partitions.txt

python3 ~/scratch/github/genbank/partitions.py -i 1_aa_partitions.txt -t aa -o raxml
python3 ~/scratch/github/genbank/partitions.py -i 2_nt_partitions.txt -t nt -o raxml

echo 'Check supermatrix and remove and taxa with no COI sequence'

# Change partitions to match a previous run (harder with codon model. Postpone!)

#sed -i -E "s/aa\///g" 1_aa_partitions.txt
#sed -i -E "s/nt\///g" 2_nt_partitions.txt
#sed -i -E "s/.fasta//g" *_partitions.txt
#sed -i -E "s/.fasta//g" ../iqtree/*.best_scheme.nex

#~/scratch/github/phylostuff/partitioner.py -p ../iqtree/nttree_MF.best_scheme.nex < 2_nt_partitions.txt > partitions_gene.nex
#~/scratch/github/phylostuff/partitioner.py -p ../iqtree/aa_MF.best_scheme.nex < 1_aa_partitions.txt > partitions_aa.nex