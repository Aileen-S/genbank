
Python
perl
MAFFT
IQTREE
PTP
RAxML

On desktop
Aliview
FigTree


catfasta2phyml.pl

get_concat_recs.py
partitions.py
ptp_get_metadata.py
ptp_filter_output.py

translate.py
backtranslate.py

### Initial GenBank Search ###

Or start with COX1 alignment/BLAST search?

python3 get_concat_recs.py -e <email> -t <taxonomy> -i both


### Constraint Tree ###

# Look at metadata and choose records for constraint tree, then write tree file.
# Look at record count for each gene. Consider excluding genes with less than a chosen number of records.


### Alignment ###

# Move metadata

mkdir metadata
mv metadata.csv metadata

# Find and replace special characters in names

for file in *
do
  sed -i 's/ /_/g' $file
  sed -i 's/\.//g' $file
  sed -i 's/\///g' $file
  sed -i 's/<//g' $file
  sed -i 's/>/ /2' $file
  sed -i 's/([^]]*)//' $file
done

# If running on desktop, try this if the above does not work
for file in *
do
  sed -i '' 's/ /_/g' $file
  sed -i '' 's/\.//g' $file
  sed -i '' 's/\///g' $file
  sed -i '' 's/<//g' $file
  sed -i '' 's/>/ /2' $file
  sed -i '' 's/([^]]*)//' $file
done

# Organise files

mkdir n1_raw n2_aa n3_aaal n4_ntal m1_raw m2_aa m3_aaal m4_ntal r1_raw r4_ntal
mv 12* 16* 18* 28* r1_raw
mv *fasta n1_raw
cd n1_raw
mv COX* CY* ATP* ND* ../m1_raw
cd ..
ls *

# Translation

for file in n1_raw/*
do
   ~/scratch/github/biotools/translate.py 1 < $file > n2_aa/${file#*/}
done

for file in m1_raw/*
do
   ~/scratch/github/biotools/translate.py 5 < $file > m2_aa/${file#*/}
done

# Alignment
# Change --globalpair to --localpair for increased accurancy if less than 200 sequences.

for file in n2_aa/*
do
mafft --6merpair --maxiterate 1000 $file >  n3_aaal/${file#*/}
done

for file in m2_aa/*
do
mafft --6merpair --maxiterate 1000 $file >  m3_aaal/${file#*/}
done

for file in r1_raw/*
do
mafft --6merpair --maxiterate 1000 $file >  r4_ntal/${file#*/}
done

# Backtranslate

for file in n3_aaal/*
do
   ~/scratch/github/biotools/backtranslate.py -s -i $file n1_raw/${file#*/} 1 > n4_ntal/${file#*/}
done

for file in m3_aaal/*
do
   ~/scratch/github/biotools/backtranslate.py -s -i $file m1_raw/${file#*/} 5 > m4_ntal/${file#*/}
done

# Backup files

mkdir phylogeny
mkdir phylogeny/aa phylogeny/nt
cp m3_aaal/* n3_aaal/* r4_ntal/* phylogeny/aa
cp m4_ntal/* n4_ntal/* r4_ntal/* phylogeny/nt

# Strip frame tags

for file in phylogeny/aa/*
do
   sed -i -E "s/;frame=[0-9]*(;$)?//" $file
done



### Copy /phylogeny/ to desktop and curate alignment files ###

# Remove unaligned sequences
# Split COX1.fasta into two halves and save as COX1a.fasta and COX1b.fasta
# Trim alignment ends where necessary


# Save GenBank accessions to gbids/gbids.txt

for file in nt/*
do
python3 ~/OneDrive\ -\ Imperial\ College\ London/github/genbank/fasta_get_gbids.py -i $file -r gbid
done

mkdir gbids
mv nt/*out gbids
cat gbids/* > gbids/all.txt
sort -u gbids/all.txt > gbids/gbids.txt
rm gbids/all.txt

# Remove GenBank accessions from fastas

for file in aa/*
do
mv $file ${file#*/}.gbids
python3 ~/scratch/github/genbank/fasta_remove_gbids.py -i ${file#*/}.gbids -o aa/${file#*/}
rm *.gbids
done

for file in nt/*
do
mv $file ${file#*/}.gbids
python3 ~/scratch/github/genbank/fasta_remove_gbids.py -i ${file#*/}.gbids -o nt/${file#*/}
rm *.gbids
done


# Build supermatrix

catfasta2phyml/catfasta2phyml.pl -c -fasta aa/* > 1_aa_supermatrix.fasta 2> 1_aa_partitions.txt
catfasta2phyml/catfasta2phyml.pl -c -fasta nt/* > 2_nt_supermatrix.fasta 2> 2_nt_partitions.txt

# Open supermatrix fastas
# Remove any record that has no sequence for COX1b (barcode region)
# Copy phylogeny/ back to server


### Phylogeny ###

# Format partition file for IQTREE

python3 ~/scratch/github/genbank/partitions.py -i 1_aa_partitions.txt -t aa -o nexus
python3 ~/scratch/github/genbank/partitions.py -i 2_nt_partitions.txt -t nt -o nexus

# Run IQTREE modelfinder and phylogeny

# Amino acid alignment with gene partitions
iqtree -s 1_aa_supermatrix.fasta -m MFP+MERGE --prefix aa_<date> -g <constraint_tree> -sp partitions_aa.txt

# Nucleotide alignment with gene partitions
iqtree -s 2_nt_supermatrix.fasta -m MFP+MERGE --prefix gene_<date> -g <constraint_tree> -sp partitions_gene.txt

# Nucleotide alignment with codon partitions
iqtree -s 2_nt_supermatrix.fasta -m MFP+MERGE --prefix codon_<date> -g <constraint_tree> -sp partitions_codon123.txt


### Species Delimitation ###

# Search for duplicate species using PTP
bPTP.py <tree_file>

# Get metadata from accessions
# Output written to metadataptp.csv
python3 ptp_get_metadata.py -e mixedupvoyage@gmail.com -f gbids.txt

# Get species list lines from output file
grep "^ " bPTP_221219.PTPhSupportPartition.txt > output.txt

# Filter PTP output lists using metadata.csv and output.txt
# Keeps all named species, and the one with most genes from each PTP grouping
# Output written to chosen.txt
python3 ptp_filter_output.py

# Refine supermatrices from PTP list
python3 ~/OneDrive\ -\ Imperial\ College\ London/github/genbank/filter_fasta.py -i 1_aa_supermatrix.fasta -f chosen.txt -o 1_aa_supermatrix.fasta
python3 ~/OneDrive\ -\ Imperial\ College\ London/github/genbank/filter_fasta.py -i 2_nt_supermatrix.fasta -f chosen.txt -o 2_nt_supermatrix.fasta


### Phylogeny ###

# Copy refined supermatrices to server
# Use IQTREE modelfinder scheme from prior step

# Amino acid alignment with gene partitions
iqtree -s 1_aa_supermatrix.fasta --prefix aa_<date> -g <constraint_tree> -Q aa_MF.best_scheme.nex

# Nucleotide alignment with gene partitions
iqtree -s 2_nt_supermatrix.fasta --prefix gene_<date> -g <constraint_tree> -Q gene_MF.best_scheme.nex

# Nucleotide alignment with codon partitions
iqtree -s 2_nt_supermatrix.fasta --prefix codon_<date> -g <constraint_tree> -Q codon_MF.best_scheme.nex

# Check trees for any unexpected results and remove tips if necessary

# Convert
