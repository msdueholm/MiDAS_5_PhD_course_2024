#!/usr/bin/bash -l
#SBATCH --job-name=MiDAS5_V13_amplicons
#SBATCH --output=job_%j_%x.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=high-mem
#SBATCH --cpus-per-task=128
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=md@bio.aau.dk

# Exit on first error and if any variables are unset
set -eu
max_threads="$(nproc)"

# Retrive and relabel V13 forward reads based in a sample list (V13_samples.txt) and merge them into one file.

mkdir V13_raw_data

while read SAMPLES
do
a="_";
NAME=$SAMPLES;
find Sequence_data/ -name $NAME$a*R1* -exec gzip -cd {} \; > V13_raw_data/$NAME.R1.fq 
usearch11 -fastx_relabel V13_raw_data/$NAME.R1.fq -prefix $NAME$a -fastqout V13_raw_data/temp.$NAME.R1.fq
cat V13_raw_data/temp.$NAME.R1.fq >> V13_raw_data/V13_forward.fastq
done < samples.txt
date

rm V13_raw_data/temp*.fq

# Usearch pipeline

mkdir V13_Usearch_analysis

### Quality filtering
usearch11 -filter_phix V13_raw_data/V13_forward.fastq -output V13_Usearch_analysis/V13_forward_phixfiltered.fq -threads $max_threads
usearch11 -fastx_truncate V13_Usearch_analysis/V13_forward_phixfiltered.fq -trunclen 250 -fastqout V13_Usearch_analysis/V13_forward_250bp.fq -threads $max_threads
usearch11 -fastq_filter V13_Usearch_analysis/V13_forward_250bp.fq -fastq_maxee 1.0 -fastaout V13_Usearch_analysis/V13_forward_250bp_filtered.fa -threads $max_threads

### Identify unique sequences
usearch11 -fastx_uniques V13_Usearch_analysis/V13_forward_250bp_filtered.fa -fastaout V13_Usearch_analysis/V13_unique.fa -sizeout -relabel uniq -threads $max_threads

### Create ASVs
usearch11 -unoise3 V13_Usearch_analysis/V13_unique.fa -zotus V13_Usearch_analysis/V13_ASV.fa -threads $max_threads
sed -i 's/Zotu/ASV/g' V13_Usearch_analysis/V13_ASV.fa 

### Create ASV table
usearch11 -otutab V13_Usearch_analysis/V13_forward_250bp_filtered.fa -zotus V13_Usearch_analysis/V13_ASV.fa -sample_delim _ -strand plus -otutabout V13_Usearch_analysis/V13_ASVtab.txt -threads $max_threads

### Predict taxonomy using sintax (MiDAS5.0_20211221)
usearch11 -sintax V13_Usearch_analysis/V13_ASV.fa -db MiDAS_5.3_sintax.fa -tabbedout V13_Usearch_analysis/V13_ASV_MiDAS_5.3.sintax -strand both -sintax_cutoff 0.8 -threads $max_threads
