#!/bin/bash
#SBATCH --job-name=blast_test
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aileen.scott@nhm.ac.uk
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-4
# Software:
# Python
# Biopython
# BLAST
# MAFFT
# CIAlign
# FastTree
# mPTP
# IQ-TREE
# RAxML NG

# Setup:
# Save taxa and NCBI taxon IDs in taxa.config config file
# Save search profile as profile.fasta
# Make sub-directory named for each taxon from config file
# Save outgroup fastas and metadata for each taxon in $taxon/outgroup

# Get taxon names from config file
taxon=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' taxa.config)
gb=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' taxa.config)

mkdir $taxon
cd $taxon


cat << p
---------------
Lab mitogenomes
---------------
p

# Remember to add outgroup
mkdir
cd lab

#python3 ~/github/genbank/get_lab_gbs.py -i ids.txt -o $taxon.gb -v gbmaster_2023-10-29

python3 ~/github/genbank/extract_lab_genes.py -g $taxon.gb -i both

mkdir frame
mv *fasta frame
mkdir noframe
mv *rf noframe

# Find reading frame
if [[ -d noframe && ! -z "$(ls -A noframe)" ]]; then
  for file in noframe/*
  do
    base=$(basename ${file%.rf})
    echo $base
    ~/scratch/github/biotools/findframe.py -r ~/voglerlab/profiles/0_AA_profiles/$base.fasta -t 5 < $file > noframe/$base.fasta
    cat  noframe/$base.fasta frame/$base.fasta > frame/$base.fasta.a
    mv frame/$base.fasta.a frame/$base.fasta
  done
fi


cat << p

------------
Add Outgroup
------------
p

# Get outgroup, if not already added with lab sequences
# IDs, 2 per family, saved in /mbl/share/workspaces/groups/voglerlab/ailes/outgroup
mkdir outgroup
cd outgroup

pushd ~/voglerlab/outgroup/
cat chosen_families > target/outgroup.gb
popd



# Buprestidae 
cat Dascillidae.gb Rhipiceridae.gb Byrrhidae.gb Dryopidae.gb Lampyridae.gb Elateridae.gb Eucinetidae.gb Clambidae.gb Derodontidae.gb Scirtidae.gb Carabidae.gb > ~/scratch/240802pipe/Buprestidae/outgroup/outgroup.gb

# Cicindellidae
cat Carabidae.gb Dytiscidae.gb Noteridae.gb Haliplidae.gb Gyrinidae.gb Trachypachidae.gb Cupedidae.gb Torridincolidae.gb > ~/scratch/240802pipe/Cicindelidae/outgroup/outgroup.gb

# Curculionidae
cat Brentidae.gb Attelabidae.gb Belidae.gb Anthribidae.gb Nemonychidae.gb Chrysomelidae.gb Cerambycidae.gb Megalopodidae.gb Erotylidae.gb Helotidae.gb Cucujidae.gb > ~/scratch/240802pipe/Curculionidae/outgroup.gb

# Elateridae
cat  Lampyridae.gb Rhagophthalmidae.gb Phengodidae.gb Lycidae.gb Cantharidae.gb Eucnemidae.gb Throscidae.gb Artematopodidae.gb Omethidae.gb Byrrhidae.gb Dryopidae.gb Buprestidae.gb  > ~/scratch/240802pipe/Elateridae/outgroup/outgroup.gb

# Lampyridae
cat Rhagophthalmidae.gb Phengodidae.gb Elateridae.gb Lycidae.gb Cantharidae.gb Eucnemidae.gb Throscidae.gb Artematopodidae.gb Omethidae.gb Byrrhidae.gb Dryopidae.gb Buprestidae.gb  > ~/scratch/240802pipe/Lampyridae/outgroup/outgroup.gb

# Tenebrionidae
cat Ciidae.gb Melandryidae.gb Zopheridae.gb Archeocrypticidae.gb Mycetophagidae.gb Meloidae.gb Anthicidae.gb Pyrochroidae.gb Salpingidae.gb Scraptiidae.gb Boridae.gb Trictenotomidae.gb Oedemeridae.gb Mycteridae.gb Aderidae.gb Mordellidae.gb Ripiphoridae.gb Lymexylidae.gb Cleridae.gb Trogossitidae.gb Melyridae.gb > ~/scratch/240802pipe/Tenebrionidae/outgroup/outgroup.gb

python3 ~/github/genbank/extract_lab_genes.py -g outgroup.gb
awk '/^LOCUS/ {print $2}' outgroup.gb > outgroup.ids

mkdir frame
mv *fasta frame
mkdir noframeawk '/^LOCUS/ {print $2}' outgroup.gb > outgroup.ids
mv *rf noframe

# Find reading frame
if [[ -d noframe && ! -z "$(ls -A noframe)" ]]; then
  for file in noframe/*
  do
    base=$(basename ${file%.rf})
    echo $base
    ~/scratch/github/biotools/findframe.py -r ~/voglerlab/profiles/0_AA_profiles/$base.fasta -t 5 < $file > noframe/$base.fasta
    cat  noframe/$base.fasta frame/$base.fasta > frame/$base.fasta.a
    mv frame/$base.fasta.a frame/$base.fasta
  done
fi

cat << p

---------------------
Search BLAST Database
---------------------
p

# Direct to BLAST database
export BLASTDB=/mnt/shared/apps/databases/ncbi/
export NCBI_API_KEY=06be85de839f1f93f5b019d354d80a9c7e09

# Get NCBI taxonomy ID for chosen taxon with script from BLAST package
output=$(get_species_taxids.sh -n $taxon)
txid=$(echo "$output" | grep 'Taxid' | awk '{print $3}')
echo "Taxon ID: $txid"

# Get species TXID list
echo "Getting txid list"
get_species_taxids.sh -t $txid > $taxon.txids

# Search BLAST nucleotide database using profile, and limit search using taxon ID list
# Output list of accessions

# Search BLAST nucleotide database using profile, limit search using taxon ID list, output list of accessions
blastn -db nt -query ~/voglerlab/profiles/blast/COX1a.fasta -taxidlist $taxon.txids -out COX1a.blast -max_target_seqs 100000 -outfmt '6 sacc'
blastn -db nt -query ~/voglerlab/profiles/blast/COX1b.fasta -taxidlist $taxon.txids -out COX1b.blast -max_target_seqs 100000 -outfmt '6 sacc'


cat << p

---------------------
Get GenBank Sequences
---------------------
p

# # Genes available with this script are ATP6, ATP8, COX1, COX2, COX3, CYTB, ND1, ND2, ND3, ND4, ND4L, ND5, ND6,
# # AK, CAD, EF1A, H3, RNApol, Wg, 12S, 16S, 18S, 28S
# python3 ~/scratch/github/genbank/search_gb_file.py -g ~/voglerlab/genbank/$gb.gb -a COX1a.blast -i both -l -c
# mv COX1.fasta COX1a.fasta
# mv metadata.csv meta_COX1a.csv
# python3 ~/scratch/github/genbank/search_gb_file.py -g ~/voglerlab/genbank/$gb.gb -a COX1b.blast -i both -l -c
# mv COX1.fasta COX1b.fasta
# mv metadata.csv meta_COX1b.csv
# python3 ~/scratch/github/genbank/search_gb_file.py -g ~/voglerlab/genbank/$gb.gb -i both -l
# rm COX1.fasta
# mv metadata.csv meta_full.csv

python3 ~/scratch/github/voglerlab/get_genbanks.py -f COX1a.blast -r gbid -e mixedupvoyage@gmail.com -c
mv COX1.fasta COX1a.fasta
mv metadata.csv meta_COX1a.csv
python3 ~/scratch/github/voglerlab/get_genbanks.py -f COX1b.blast -r gbid -e mixedupvoyage@gmail.com -c
mv COX1.fasta COX1b.fasta
mv metadata.csv meta_COX1b.csv
python3 ~/scratch/github/voglerlab/get_genbanks.py -t $gb -e mixedupvoyage@gmail.com
rm COX1.fasta
mv metadata.csv meta_full.csv

# Combine metadata
cat meta* > meta1.csv
sort meta1.csv | uniq > meta2.csv
awk '/^ncbi/{print; next} {a[NR]=$0} END{for(i in a) print a[i]}' meta2.csv > metadata.csv
rm meta1.csv meta2.csv

# Get frame tags if necessary
mkdir frame
mv *fasta frame
mkdir noframe
mv *rf noframe

# Find reading frame
if [[ -d noframe && ! -z "$(ls -A noframe)" ]]; then
  for file in noframe/*
  do
    base=$(basename $file)
    echo $base
    ~/scratch/github/biotools/findframe.py -r ~/voglerlab/profiles/0_AA_profiles/${base#*/} -t 5 < $file > noframe/${base%rf}fasta
    cat $file frame/$base > frame/$base.a
    mv frame/$base.a frame/$base
  done
fi



cat << p

-----------
BOLD search
-----------
p

mkdir bold
cd bold

cp ~/scratch/240611tiger/$taxon/bold/raw_meta* .

# Search BOLD database online (or use the same script to search a previously downloaded BOLD file
#Rscript ~/scratch/github/voglerlab/bold_cli.R -t $taxon -m ../genbank/metadata.csv

Rscript ~/scratch/github/voglerlab/bold_cli.R -c raw_metadata.tsv -m ../genbank/metadata.csv -g


cat << p

-----------------------
Get frame tags for BOLD
-----------------------
p

mkdir noframe
mkdir frame
mv *fasta noframe

# Find reading frame
for file in noframe*
do
  base=$(basename $file)
  echo $base
  findframe.py -r ~/voglerlab/profiles/0_AA_profiles/$base -t 5 < $file > frame/$base
done

# Copy RNA fastas to 'frame' directory
cp noframe/1* noframe/2* frame
cd ..



cat << p

-------------------
Merge raw sequences
-------------------
p

mkdir alignment
cd alignment

# Define the list of possible file names
file_names=("12S.fasta" "16S.fasta" "18S.fasta" "28S.fasta" "AK.fasta" "ATP6.fasta" "ATP8.fasta" "CAD.fasta" "COX1.fasta" "COX1a.fasta" "COX1b.fasta" "COX2.fasta" "COX3.fasta" "CYTB.fasta" "EF1A.fasta" "H3.fasta" "ND1.fasta" "ND2.fasta" "ND3.fasta" "ND4.fasta" "ND4L.fasta" "ND5.fasta" "ND6.fasta" "Wg.fasta")

# List of directories to check
directories=(
  "../lab/frame"
  "../genbank/frame"
  "../bold/frame"
  "../outgroup/frame"
)

# Process each file name
for f in "${file_names[@]}"; do
  output="$f"
  : > "$output"  # Create or truncate the output file

  for dir in "${directories[@]}"; do
    if [ -e "$dir/$f" ]; then
      cat "$dir/$f" >> "$output"
    fi
  done
done

# Split into RNA, mitogenome protein coding and nuclear protein coding
mkdir n1_raw m1_raw r1_raw
mkdir n2_aa n3_aaal n4_ntal m2_aa m3_aaal m4_ntal r4_ntal
mv 12S.fasta 16S.fasta 18S.fasta 28S.fasta r1_raw
mv COX*fasta CYTB.fasta ATP*fasta ND*fasta m1_raw
mv RNApol.fasta Wg.fasta EF1A.fasta H3.fasta CAD.fasta AK.fasta n1_raw

cat << p

--------------
COX1 sequences
--------------
p
# Merge COX1a and COX1b, removing duplicated. Seperate manually after alignment
# Will not work if some sequences have already been separated and have the same IDs
mkdir COI
mv m1_raw/COX1* COI
cat COI/COX1* > COI/COX1all.fasta
~/scratch/github/voglerlab/fasta_dups_length.py -i COI/COX1all.fasta -o m1_raw/COX1.fasta -d


cat << p

--------------------
Translate to Protein
--------------------
p

# Mitochondrial genes
for file in m1_raw/*
do
  echo 'translate' $file
   ~/scratch/github/biotools/translate.py 5 < $file > m2_aa/${file#*/}
done

# Nuclear genes
for file in n1_raw/*
do
  echo 'translate' $file
   ~/scratch/github/biotools/translate.py 1 < $file > n2_aa/${file#*/}
done


cat << p

----------------
Align to Profile
----------------
p

# Mitochondrial
for file in m2_aa/*
do
  echo $file
  mafft --add $file --maxiterate 2 --adjustdirection --thread 12 ~/voglerlab/profiles/0_AA_profiles/${file#*/} >  m3_aaal/${file#*/}
done

# Nuclear
for file in n2_aa/*
do
  echo $file
  mafft --add $file --maxiterate 2 --adjustdirection --thread 12 ~/voglerlab/profiles/0_AA_profiles/${file#*/} >  n3_aaal/${file#*/}
done

# RNA
for file in r1_raw/*
do
  echo $file
  mafft --add $file --maxiterate 2 --adjustdirection --thread 12 ~/voglerlab/profiles/0_AA_profiles/${file#*/} >  r4_ntal/${file#*/}
done


cat << p

--------------
Remove Profile
--------------
p

cd m3_aaal
for file in *
do
  echo $file
  mv $file ${file}_bak
  perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' ${file}_bak | grep -A 1 -P '^>(?!PROFILE::).*' > $file
  rm ${file}_bak
done
cd ..

cd n3_aaal
for file in *
do
  echo $file
  mv $file ${file}_bak
  perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' ${file}_bak | grep -A 1 -P '^>(?!PROFILE::).*' > $file
  rm ${file}_bak
done
cd ..

cd r4_ntal
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

# Mitochondrial
for file in m3_aaal/*
do
  echo $file
  ~/scratch/github/biotools/backtranslate.py -i $file m1_raw/${file#*/} 5 > m4_ntal/${file#*/}
done

# Nuclear
for file in n3_aaal/*
do
  echo $file
  ~/scratch/github/biotools/backtranslate.py -i $file n1_raw/${file#*/} 5 > n4_ntal/${file#*/}
done

# Tidy up sequence ends and remove duplicates
for file in m4_ntal/*
do
  base=$(basename $file)
  ~/scratch/github/voglerlab/fasta_dups_length.py -i $file -o m4_ntal/$base.f -l -d
  mv m4_ntal/$base.f m4_ntal/$base
done

for file in n4_ntal/*
do
  base=$(basename $file)
  ~/scratch/github/voglerlab/fasta_dups_length.py -i $file -o n4_ntal/$base.f -l -d
  mv n4_ntal/$base.f n4_ntal/$base
done

for file in r4_ntal/*
do
  base=$(basename $file)
  ~/scratch/github/voglerlab/fasta_dups_length.py -i $file -o r4_ntal/$base.f -l -d
  mv r4_ntal/$base.f r4_ntal/$base
done


cat << p

----------------
Clean Alignments
----------------
p

mkdir m5_aa_cia m6_nt_cia n5_aa_cia n6_nt_cia r6_nt_cia

for file in m3_aaal/*; do
  echo $file
  base=$(basename $file)
  if [[ $base == ATP8.fasta ]]; then
    CIAlign --infile $file --outfile_stem m5_aa_cia/${base%.fasta} --remove_divergent --remove_divergent_minperc 0.5
  else
    CIAlign --infile $file --outfile_stem m5_aa_cia/${base%.fasta} --remove_divergent --remove_divergent_minperc 0.5 --crop_divergent --crop_divergent_min_prop_nongap 0.1 --remove_insertions
  fi
done

for file in m4_ntal/*; do
  echo $file
  base=$(basename $file)
  if [[ $base == ATP8.fasta ]]; then
    CIAlign --infile $file --outfile_stem m6_nt_cia/${base%.fasta} --remove_divergent --remove_divergent_minperc 0.6 --crop_divergent --crop_divergent_min_prop_nongap 0.1 --remove_insertions
  else
    CIAlign --infile $file --outfile_stem m6_nt_cia/${base%.fasta} --remove_divergent --remove_divergent_minperc 0.65 --crop_divergent --crop_divergent_min_prop_nongap 0.1 --remove_insertions
  fi
done
 
for file in n3_aaal/*
do
  echo $file
  base=$(basename $file)
  CIAlign --infile $file --outfile_stem n5_aa_cia/${base%.fasta} --remove_divergent --remove_divergent_minperc 0.5 --crop_divergent --crop_divergent_min_prop_nongap 0.1 --remove_insertions
done

for file in n4_ntal/*
do
  echo $file
  base=$(basename $file)
  CIAlign --infile $file --outfile_stem n6_nt_cia/${base%.fasta} --remove_divergent --remove_divergent_minperc 0.65 --crop_divergent --crop_divergent_min_prop_nongap 0.1 --remove_insertions
done

for file in r4_ntal/*
do
  echo $file
  base=$(basename $file)
  CIAlign --infile $file --outfile_stem r6_nt_cia/${base%.fasta} --remove_divergent --remove_divergent_minperc 0.65 --crop_divergent --crop_divergent_min_prop_nongap 0.1 --remove_insertions
done

# Split COX1

cd m6_nt_cia
awk '
/^>/ {  # When encountering a header line
    if (NR > 1) {  # For all except the first header line
        if (length(seq) > 0 && !is_all_gaps(seq)) {  # Only process non-empty sequences that are not all gaps
            print_seq("COX1a.fasta", "COX1b.fasta", header, seq)
        }
    }
    header = $0  # Save the header
    seq = ""  # Reset the sequence
    next
}
{
    seq = seq $0  # Append sequence lines
}
END {
    if (length(seq) > 0 && !is_all_gaps(seq)) {  # Process the last sequence if it is not all gaps
        print_seq("COX1a.fasta", "COX1b.fasta", header, seq)
    }
}

function print_seq(file1, file2, header, seq,    len1, len2) {
    len1 = 720
    len2 = length(seq)
    if (len2 > len1) {
        if (length(substr(seq, 1, len1)) > 0 && !is_all_gaps(substr(seq, 1, len1))) {  # Only write non-empty sequences that are not all gaps
            print header > file1
            print substr(seq, 1, len1) > file1
        }
        if (length(substr(seq, len1 + 1)) > 0 && !is_all_gaps(substr(seq, len1 + 1))) {  # Only write non-empty sequences that are not all gaps
            print header > file2
            print substr(seq, len1 + 1) > file2
        }
    } else {
        if (length(seq) > 0 && !is_all_gaps(seq)) {  # Only write non-empty sequences that are not all gaps
            print header > file1
            print seq > file1
        }
    }
}

function is_all_gaps(seq,    i) {
    for (i = 1; i <= length(seq); i++) {
        if (substr(seq, i, 1) != "-") {
            return 0  # Not all gaps
        }
    }
    return 1  # All gaps
}
' COX1_cleaned.fasta

cd ../m5_aa_cia
awk '
/^>/ {  # When encountering a header line
    if (NR > 1) {  # For all except the first header line
        if (length(seq) > 0 && !is_all_gaps(seq)) {  # Only process non-empty sequences that are not all gaps
            print_seq("COX1a.fasta", "COX1b.fasta", header, seq)
        }
    }
    header = $0  # Save the header
    seq = ""  # Reset the sequence
    next
}
{
    seq = seq $0  # Append sequence lines
}
END {
    if (length(seq) > 0 && !is_all_gaps(seq)) {  # Process the last sequence if it is not all gaps
        print_seq("COX1a.fasta", "COX1b.fasta", header, seq)
    }
}

function print_seq(file1, file2, header, seq,    len1, len2) {
    len1 = 250
    len2 = length(seq)
    if (len2 > len1) {
        if (length(substr(seq, 1, len1)) > 0 && !is_all_gaps(substr(seq, 1, len1))) {  # Only write non-empty sequences that are not all gaps
            print header > file1
            print substr(seq, 1, len1) > file1
        }
        if (length(substr(seq, len1 + 1)) > 0 && !is_all_gaps(substr(seq, len1 + 1))) {  # Only write non-empty sequences that are not all gaps
            print header > file2
            print substr(seq, len1 + 1) > file2
        }
    } else {
        if (length(seq) > 0 && !is_all_gaps(seq)) {  # Only write non-empty sequences that are not all gaps
            print header > file1
            print seq > file1
        }
    }
}

function is_all_gaps(seq,    i) {
    for (i = 1; i <= length(seq); i++) {
        if (substr(seq, i, 1) != "-") {
            return 0  # Not all gaps
        }
    }
    return 1  # All gaps
}
' COX1_cleaned.fasta
cd ..


cat << p

---------------------------------
Prepare sequences for supermatrix
---------------------------------
p

mkdir ../supermatrix
mkdir ../supermatrix/aa
mkdir ../supermatrix/nt

cp m5_aa_cia/*fasta n5_aa_cia/*fasta r6_nt_cia/*fasta ../supermatrix/aa
cp m6_nt_cia/*fasta n6_nt_cia/*fasta r6_nt_cia/*fasta ../supermatrix/nt

cd ../supermatrix

# Delete files with less than 10 sequences

for file in aa/*
do
  sequence_count=$(grep -c "^>" "$file")
  # Check if the sequence count is less than 10
  if (( sequence_count < 10 )); then
    echo "Deleting $file, as it contains only $sequence_count sequences"
    rm "$file"
  fi
done

for file in nt/*
do
  sequence_count=$(grep -c "^>" "$file")
  # Check if the sequence count is less than 10
  if (( sequence_count < 10 )); then
    echo "Deleting $file, as it contains only $sequence_count sequences"
    rm "$file"
  fi
done

# Remove frame tags
for file in aa/*
do
   sed -i -E "s/;frame=[0-9]*(;$)?//" $file
done
echo "Frame tags removed from files in aa"

for file in nt/*
do
   sed -i -E "s/;frame=[0-9]*(;$)?//" $file
done
echo "Frame tags removed from files in nt"

# Remove _R_ from reversed sequences
for file in aa/*
do
   sed -i -E "s/^>_R_/>/" $file
done

for file in nt/*
do
   sed -i -E "s/^>_R_/>/" $file
done


cat << p

----------------
Combine metadata
----------------
p

cd ..

Rscript ~/github/test/combine_metadata.R -b bold/metadata.csv -g genbank/metadata.csv -l lab/ids -o outgroup/ids -m ~/scratch/SITE-100_Database_mitogenomes_240731.csv -c supermatrix/rename.csv

cd supermatrix

# Rename sequences with taxonomy
mkdir ids_aa ids_nt
for file in nt/*
do
  base=$(basename $file)
  python3 ~/scratch/github/genbank/rename_fasta.py -i $file -c rename.csv -r ids_nt/${base%fasta}csv -o nt/$base.r
  mv nt/${base#*/}.r nt/${base#*/}
done

for file in aa/*
do
  base=$(basename $file)
  python3 ~/scratch/github/genbank/rename_fasta.py -i $file -c rename.csv -r ids_aa/${base%fasta}csv -o aa/$base.r
  mv aa/${base#*/}.r aa/${base#*/}
done

cat << p

----------------------------
Remove sequences without COI
----------------------------
p

mkdir aa_coi nt_coi

for file in aa/*
do
  echo $file
  python3 ~/scratch/github/genbank/filter_fasta.py -i $file -s aa/COX1_cleaned.fasta -t fasta -f aa_coi/${file#*/}
  # Remove gap only columns
  echo 'Removing gap-only columns'
  ~/scratch/github/fastagap/degap_fasta_alignment.pl aa_coi/${file#*/} > aa_coi/${file#*/}.g
  mv aa_coi/${file#*/}.g aa_coi/${file#*/}
done

for file in nt/*
do
  echo $file
  python3 ~/scratch/github/genbank/filter_fasta.py -i $file -s aa/COX1_cleaned.fasta -t fasta -f nt_coi/${file#*/}
  # Remove gap only columns
  echo 'Removing gap-only columns'
  ~/scratch/github/fastagap/degap_fasta_alignment.pl nt_coi/${file#*/} > nt_coi/${file#*/}.g
  mv nt_coi/${file#*/}.g nt_coi/${file#*/}
done

rm */COX1_cleaned.fasta


# Remove empty files
find aa_coi -type f -empty -delete
find nt_coi -type f -empty -delete


cat << p

---------------------------------------
Build Supermatrices and Partition Files
---------------------------------------
p

~/scratch/github/catfasta2phyml/catfasta2phyml.pl -c -fasta aa_coi/* > 1_aa_coi_supermatrix.fasta 2> 1_aa_coi_partitions.txt
~/scratch/github/catfasta2phyml/catfasta2phyml.pl -c -fasta nt_coi/* > 2_nt_coi_supermatrix.fasta 2> 2_nt_coi_partitions.txt

python3 ~/scratch/github/genbank/partitions.py -i 1_aa_coi_partitions.txt -t aa -o raxml
python3 ~/scratch/github/genbank/partitions.py -i 2_nt_coi_partitions.txt -t nt -o raxml

cd ..

cat << p

---------
Phylogeny
---------
p

mkdir mptp

# Run FastTree
raxml-ng --search1 --msa supermatrix/2_nt_coi_supermatrix.fasta --prefix mptp/search1 --model supermatrix/partitions_gene.txt --threads 8

cat << p

-------------------------
mPTP Species Delimitation
-------------------------
p

cd mptp
# Get outgroup
outgroup=$(awk 'NR==1' outgroup.txt)

# Calculate minbr value
echo "Getting minbr value"
mptp --minbr_auto ../supermatrix/2_nt_coi_supermatrix.fasta --tree_file search1.raxml.bestTree --output_file minbr  --outgroup $outgroup

# Copy minbr value from slurm output
minbr=$(awk 'END {print $NF}' minbr.txt)
echo "minbar value = " $minbr

# Run PTP
echo "Running mPTP"
mptp --ml --single --minbr $minbr --tree_file search1.raxml.bestTree --output_file mptp  --outgroup $outgroup

# Run MCMC sampling to determine statistical significance
#echo "
#Running MCMC"
#mptp --mcmc 1000000 --mcmc_sample 1000 --mcmc_log 1000 --mcmc_runs 2 --multi --minbr $minbr --tree_file ../trees/$taxon.raxml.bestTree --output_file mcmc
#cd ..

# Filter taxa
# Choose taxon with most nucleotides for each PTP species
python3 ~/scratch/github/genbank/ptp_filter_output.py -i mptp.txt -s ../supermatrix/2_nt_coi_supermatrix.fasta -o ptp_chosen.txt
# Filter supermatrix
python3 ~/scratch/github/genbank/filter_fasta.py -i ../supermatrix/2_nt_coi_supermatrix.fasta -f ../supermatrix/2_nt_coi_ptp_supermatrix.fasta -s ptp_chosen.txt -t list
python3 ~/scratch/github/genbank/filter_fasta.py -i ../supermatrix/1_aa_coi_supermatrix.fasta -f ../supermatrix/1_aa_coi_ptp_supermatrix.fasta -s ptp_chosen.txt -t list

cd ..

cat << p

-----
RAxML
-----
p

mkdir raxml
cd raxml

# Get random number seed for RAxML and save to slurm output
seed=$RANDOM
echo 'RAxML random number seed = '$seed

# Run RaxML
raxml-ng --search1 --msa ../supermatrix/2_nt_coi_ptp_supermatrix.fasta --model ../supermatrix/partitions_gene.txt --prefix nt_gene --threads 8 --seed $seed

run_treeshrink.py -t nt_gene.raxml.bestTree -o treeshrink

# Filter supermatrix
python3 ~/scratch/github/genbank/filter_fasta.py -i ../supermatrix/2_nt_coi_ptp_supermatrix.fasta -f ../supermatrix/2_nt_coi_ptp_ts_supermatrix.fasta -s treeshrink/output.bestTree -t tree
python3 ~/scratch/github/genbank/filter_fasta.py -i ../supermatrix/1_aa_coi_ptp_supermatrix.fasta -f ../supermatrix/1_aa_coi_ptp_ts_supermatrix.fasta -s treeshrink/output.bestTree -t tree
python3 ~/scratch/github/genbank/ry_code.py -i ../supermatrix/2_nt_coi_ptp_supermatrix.fasta -o ../supermatrix/3_ry_coi_ptp_supermatrix.fasta


cat << p

--------
RAxML ML
--------
p

#!/bin/bash
#SBATCH --mem=5G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-20
#SBATCH -p long

seed=$RANDOM
echo "Random Seed: $RANDOM

mkdir ml_nt_gene ml_codon ml_aa ml_ry_gene
raxml-ng --search1 --msa ../supermatrix/2_nt_coi_ptp_supermatrix.fasta --model ../supermatrix/partitions_gene.txt --prefix ml_nt_gene/$SLURM_ARRAY_TASK_ID --threads 8 --seed $seed --force perf_threads
raxml-ng --search1 --msa ../supermatrix/2_nt_coi_ptp_supermatrix.fasta --model ../supermatrix/partitions_codon123.txt --prefix ml_nt_codon/$SLURM_ARRAY_TASK_ID --threads 8 --seed $seed --force perf_threads

raxml-ng --search1 --msa ../supermatrix/1_aa_coi_ptp_supermatrix.fasta --model ../supermatrix/partitions_aa.txt --prefix ml_aa_gene/$SLURM_ARRAY_TASK_ID --threads 8 --seed $seed --force perf_threads

#raxml-ng --search1 --msa ../supermatrix/3_ry_coi_ptp_supermatrix.fasta --model ../supermatrix/partitions_gene.txt --prefix ml_ry_gene/$SLURM_ARRAY_TASK_ID --threads 8 --seed $seed --force perf_threads
#raxml-ng --search1 --msa ../supermatrix/3_ry_coi_ptp_supermatrix.fasta --model ../supermatrix/partitions_codon123.txt --prefix ml_ry_codon/$SLURM_ARRAY_TASK_ID --threads 8 --seed $seed --force perf_threads


# Then find best trees when all finished
# Check they're finished
grep bestTree ml_nt_gene/*log | wc -l
grep bestTree ml_codon/*log | wc -l

# Find best trees
best=$(grep 'Final LogLikelihood' ml_nt_gene/*log | awk -F '[:,]' '{print $1, $3}' | sort -k2,2nr | head -n 1 | awk -F '.' '{print $1}')
cp $best.raxml.bestTree ml_nt_gene.bestTree

best=$(grep 'Final LogLikelihood' ml_codon/*log | awk -F '[:,]' '{print $1, $3}' | sort -k2,2nr | head -n 1 | awk -F '.' '{print $1}')
cp $best.raxml.bestTree ml_codon.bestTree

best=$(grep 'Final LogLikelihood' ml_aa/*log | awk -F '[:,]' '{print $1, $3}' | sort -k2,2nr | head -n 1 | awk -F '.' '{print $1}')
cp $best.raxml.bestTree ml_aa.bestTree
