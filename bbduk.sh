#!/bin/bash
#need to change this depending on the file names
for i in $(ls *.fastq.gz | rev | cut -c 16- | rev | uniq); do
    /home/nicolekirchoff/tools/bbmap/bbduk.sh -Xmx20g \
    in1=${i}R1_001.fastq.gz \
    in2=${i}R2_001.fastq.gz \
    out=/media/nicolekirchoff/DATA/Projects/Thomas_vaccine/fastqs_trimmed/${i}R1_001.fastq.gz \
    out2=/media/nicolekirchoff/DATA/Projects/Thomas_vaccine/fastqs_trimmed/${i}R2_001.fastq.gz \
    ref=/home/nicolekirchoff/tools/bbmap/resources/adapters.fa \
    ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=10 \
    stats=/media/nicolekirchoff/DATA/Projects/Thomas_vaccine/fastqs_trimmed/bbduk_stats/${i}AT_stats.txt
done
