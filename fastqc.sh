#!/bin/bash
#Script runs FastQC on all fastqs in the input directory

# Report:
echo -e "\n## Starting FastQC script."
date

# SETUP --------------------------------------------------------
# Command-line arguments:
indir=$1
outdir=$2

# Create output directory if it doesn't already exist:
mkdir -p $outdir

echo "## Input dir:       $indir"
echo "## Output dir:      $outdir"

# RUNNING FASTQC ---------------------------------------------

fastqc $indir/*fastq.gz -t 2 -o $outdir
