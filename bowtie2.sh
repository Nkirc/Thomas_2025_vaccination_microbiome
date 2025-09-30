#!/usr/bin/env bash
##Script runs a directory of fastq files through bowtie2 and return .sam files. Fastqs need to be paired end with R1 and R2 in the file name. Modified slightly from script by instructors a$

# Report:
echo -e "\n## Starting bowtie2 script."
date

# SETUP --------------------------------------------------------
# Command-line arguments:
indir=$1
outdir=$2

# Create output directory if it doesn't already exist:
mkdir -p $outdir

echo "## Input dir:       $indir"
echo "## Output dir:      $outdir"

# RUN BOWTIE2 --------------------------------------------------------

echo -e "\n\n## Running bowtie2..."

for i in $(ls *.gz | rev | cut -c 16- | rev | uniq) ; do \
        echo -e "## Running sample ${i} \n\n"
        bowtie2 -p 40 -x /home/nicolekirchoff/RefSeq/Mus_musculus_genome/ncbi_dataset/data/GCF_000001635.27/mouse \
        -1 $indir/${i}R1_001.fastq.gz -2 $indir/${i}R2_001.fastq.gz \
        -S $outdir/${i}.sam 2> "${i}stats.txt"

        echo -e "\n\n ------------------------------------------------------------\n"

done
