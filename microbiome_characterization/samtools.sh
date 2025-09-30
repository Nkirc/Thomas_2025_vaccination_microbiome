#!/usr/bin/env bash

##Script runs a directory of .sam files from bowtie2 and outputs 2 .bam files: one file with the reads that mapped to the reference genome and one file with reads the did NOT map to the reference genome.
#Hardcoded portion below needs changing based on filenames

# Report:
echo -e "\n## Starting bowtie2 script."
date

# SETUP --------------------------------------------------------
# Command-line arguments:
indir=$1
unmap_outdir=$2
map_outdir=$3

# Create output directory if it doesn't already exist:
mkdir -p $unmap_outdir
mkdir -p $map_outdir

echo "## Input dir:       $indir"
echo "## Unmapped output dir:      $unmap_outdir"
echo "## Mapped output dir:      $map_outdir"

# RUN BOWTIE2 --------------------------------------------------------

echo -e "\n\n## Running samtools..."

for i in $(ls *.sam | rev | cut -c 5- | rev | uniq); do \
        echo -e "## Running sample ${i} \n\n"
       samtools view -b -f 12 -@ 30 $indir/${i}.sam -o $unmap_outdir/${i}unmapped.bam \
       #samtools view -b -F 12 -@ 30 $indir/${i}.sam.gz -o $map_outdir/${i}mapped.bam \

        echo -e "\n\n ------------------------------------------------------------\n"

done
