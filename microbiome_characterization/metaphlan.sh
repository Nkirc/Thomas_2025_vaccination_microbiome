conda activate metaphlan

for i in *fastq.gz; do metaphlan $i --input-type fastq -nproc 40 > ${i%.fastq.gz}_profile.txt;done

merge_metaphlan_tables.py *_profile.txt > merged_abundance_table.txt
