# First get database to search, and profile to search with.

# If using slurm, get taxon names from config file
taxon=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' config.txt)

mkdir $taxon
cd $taxon

# Format profile sequences for HMMER
hmmbuild $taxon.profile.hmm $taxon.profile.fasta

# HMMER search
# Add -A $taxonout.sth for sequence file output
nhmmer --tblout $taxon.hmmhits.txt $taxon.profile.hmm /home/ascott/voglerlab/barcodes/fastas/Adephaga/$taxon.fasta

# Get list of accessions from HMMER output file
python3 ~/scratch/github/genbank/hmm_get_accessions.py -i $taxon.hmmhits.txt -o $taxon.accessionlist.txt

# Get mitogenome sequences from HMMER hits list - longest sequence for each taxon ID is saved
mkdir raw
cd raw
python3 ~/scratch/github/genbank/get_concat_recs.py -f $taxon.accessionlist.txt -r gbid -i both -e <your_email_address>
mv metadata.csv ..
cd ..

# Translate
mkdir aa
for file in raw/*
do
  ~/scratch/github/biotools/translate.py 1 < $file > aa/${file#*/}
done
cd ..

# Align
mkdir aa_align
for file in aa/*
do
mafft --add $file --global --maxiterate 1000 --anysymbol --thread 10 ~/scratch/profiles/0_AA_profiles >  aa_align/${file#*/}
done

mkdir nt_align
for file in aa_align/*
do
  ~/scratch/github/biotools/backtranslate.py -i $file raw/${file#*/} 1 > nt_align/${file#*/}
done