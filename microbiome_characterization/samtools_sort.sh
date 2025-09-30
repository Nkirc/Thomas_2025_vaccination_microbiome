#!usr/bin/env bash

#Script sorts bam files by identifier name

#Report:
echo -e "\n## Starting samtools script."

#SETUP---------------------------------------------
#Command-line arguments
indir=$1
outdir=$2

echo "## Input directory:	$indir"
echo "## Output directory:	$outdir"

echo"\n## Running samtools.."


for i in $(ls *unmapped.bam | rev | cut -c 13- | rev | uniq); do \
        echo -e "## Running sample ${i} \n\n"
	  samtools sort -@ 30 -n $indir/${i}unmapped.bam -o $outdir/${i}unmap_sort.bam
        echo -e "\n\n ------------------------------------------------------------\n"

done
