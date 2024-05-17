#!/usr/bin/bash -l
#SBATCH --job-name=MiDAS5_V4_amplicons
#SBATCH --output=job_%j_%x.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=general
#SBATCH --cpus-per-task=128
#SBATCH --mem=20G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=md@bio.aau.dk

# Exit on first error and if any variables are unset
set -eu
max_threads="$(nproc)"

# Retrive and relabel V4 forward reads based in a sample list (samples.txt) and merge them into one file.

mkdir V4_raw_data

while read SAMPLES
do
a="_";
NAME=$SAMPLES;
find Sequence_data/ -name $NAME$a*R1* -exec gzip -cd {} \; > V4_raw_data/$NAME.R1.fq 
usearch11 -fastx_relabel V4_raw_data/$NAME.R1.fq -prefix $NAME$a -fastqout V4_raw_data/temp.$NAME.R1.fq
cat V4_raw_data/temp.$NAME.R1.fq >> V4_raw_data/V4_forward.fastq
done < samples.txt
date

rm V4_raw_data/temp*.fq

# Usearch pipeline
mkdir V4_Usearch_analysis

### Adapter trimming
module load cutadapt/2.8-foss-2018a-Python-3.6.4
cutadapt -g ^GGACTACHVGGGTWTCTAAT...TTACCGCGGCKGCTGGCAC -o V4_Usearch_analysis/V4_trimmed.fq --discard-untrimmed V4_raw_data/V4_forward.fastq -j $max_threads

### Quality filtering
usearch11 -fastx_revcomp V4_Usearch_analysis/V4_trimmed.fq -fastqout V4_Usearch_analysis/V4_trimmed_rc.fq -threads $max_threads
usearch11 -fastq_filter V4_Usearch_analysis/V4_trimmed_rc.fq -fastq_maxee 1.0 -fastaout V4_Usearch_analysis/V4_qcfiltered.fa -threads $max_threads

### Identify unique sequences
usearch11 -fastx_uniques V4_Usearch_analysis/V4_qcfiltered.fa -fastaout V4_Usearch_analysis/V4_unique.fa -sizeout -relabel uniq

### Create ASVs
usearch11 -unoise3 V4_Usearch_analysis/V4_unique.fa -zotus V4_Usearch_analysis/V4_ASV.fa
sed -i 's/Zotu/ASV/g' V4_Usearch_analysis/V4_ASV.fa

### Create ASV table
usearch11 -otutab V4_Usearch_analysis/V4_qcfiltered.fa -zotus V4_Usearch_analysis/V4_ASV.fa -sample_delim _ -strand plus -otutabout V4_Usearch_analysis/V4_ASVtab.txt -threads $max_threads

### Predict taxonomy using sintax (SILVA 138.1 SSURef NR99 and MiDAS5.3_)
usearch11 -sintax V4_Usearch_analysis/V4_ASV.fa -db MiDAS_5.3_sintax.fa -tabbedout V4_Usearch_analysis/V4_ASV_MiDAS_5.3.sintax -strand both -sintax_cutoff 0.8 -threads $max_threads
