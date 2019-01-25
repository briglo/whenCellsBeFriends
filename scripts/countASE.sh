#!/bin/bash
#usage 1=path/to/bam 2= path/to/vcf 3=sampleName

samtools mpileup -v -u -t AD -f /share/ClusterShare/biodata/contrib/briglo/hg19/hs37d5.fa -l "$2" -x "$1" | awk 'BEGIN{FS="\t";OFS="\t"} /^[^#]/ { chrom=$1; pos=$2; ref=$4; alts=$5; datas=$10; split(alts,alt, /,/); split(datas, data_fields, /:/); split(data_fields[2], depths, /,/); print chrom, pos, ref, alt[1], depths[1], depths[2]}' > "$3"_refAltCounts.tsv