#!/usr/bin/env bash

#Script takes R1 and R2 reads files and merges them into one read file by overlapping the two reads.
# Report:
echo -e "\n## Starting flash merging script."
date

# SETUP --------------------------------------------------------
# Command-line arguments:
indir=$1
outdir=$2
# Create output directory if it doesn't already exist:
mkdir -p $outdir
echo "## Input dir:       $indir"
echo "## Output dir:      $outdir"

# MERGE WITH FLASH --------------------------------------------------------
echo -e "\n\n## Running flash..."

for i in $(ls *.fastq.gz | rev | cut -c 12- | rev | uniq) ; do \
        echo -e "## Running sample ${i} \n\n"
	
	flash $indir/${i}R1.fastq.gz $indir/${i}R2.fastq.gz -M 150 -o ${i} -d ../fastqs_unmapped_merged/

        echo -e "\n\n ------------------------------------------------------------\n"

done
