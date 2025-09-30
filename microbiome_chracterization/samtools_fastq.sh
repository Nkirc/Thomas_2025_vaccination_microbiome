#!/usr/bin/env bash

#Script converts a sorted bam file into two fastq files (one for each of the pair)

#Report:
echo -e "\n## Starting samtools script."

#SETUP---------------------------------------------
#Command-line arguments
indir=$1
outdir=$2

echo "## Input directory:       $indir"
echo "## Output directory:      $outdir"

echo"\n## Running samtools.."


for i in $(ls *.bam | rev | cut -c 15- | rev | uniq); do \
        echo -e "## Running sample ${i} \n\n"
	#-s dev/null removes any singletons
	samtools fastq -@ 8 $indir/${i}unmap_sort.bam \
    	-1 $outdir/${i}R1.fastq.gz \
    	-2 $outdir/${i}R2.fastq.gz \
    	#-0 /dev/null -s /dev/null -n
        echo -e "\n\n ------------------------------------------------------------\n"

done
