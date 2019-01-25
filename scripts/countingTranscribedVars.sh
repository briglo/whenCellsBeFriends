sed -i "1s/^.*$/##fileformat=VCFv4.1/" IPSCMedia_Genotyping_Sample1.vcf
sed -i "1s/^.*$/##fileformat=VCFv4.1/" IPSCMedia_Genotyping_Sample2.vcf

#trim VCF for exonic (dunno if this makes it better)
module load pethum/bedtools/gcc-4.4.6/2.25.0
for i in *.vcf ; do \
qsub -V -cwd -b y -j y -N bt_"$i" -pe smp 1 \
"bedtools intersect -wa -header \
-a "$i" \
-b ../annotations/hg19_ens_exon_reduce.mapping_sorted.bed > ./subset"$i"" ; done

#run demuxlet
module load briglo/demuxlet/1.9
#see https://github.com/statgen/demuxlet for help
demuxlet --sam ./data/cellranger/IPSCMedia_scRNA_Sample2_V3/outs/possorted_genome_bam.bam --vcf ./data/vcf/subsetIPSCMedia_Genotyping_Sample2.vcf --out ./demux_
demuxlet --sam ./data/cellranger/IPSCMedia_scRNA_Sample2_V3/outs/possorted_genome_bam.bam --vcf ./data/vcf/IPSCMedia_Genotyping_Sample2.vcf --out ./output/demuxlet/AllImputedSNP

##run scrublet?
module load briglo/miniconda/3
source activate scrublet




##garve up on whats below jan 2019 coz too many people have done it and i dont care. Ill just find someone who has
 for i in subset* ; do cat "$i" | awk 'BEGIN {FS="\t";OFS="\t"} ($10 ~ /^0\|1:/) {print $1,$2-1,$2,$4,$5}' > exon_het_"$i".bed ; done

for i in subset* ; do cat "$i" | awk 'BEGIN {FS="\t";OFS="\t"} ($1 ~ /#/)  {print $1,$2-1,$2,$4,$5}' > exon_all_"$i".bed ; done


 ../

# this is now called countASE.sh
#samtools mpileup -v -u -t AD -f /share/ClusterShare/biodata/contrib/briglo/hg19/hs37d5.fa -l het_subsetSYD-40370185.dedup.realigned.recalibrated.hc.gvcf.gz.bed -x briglo/RNA/6_CIRCA_MIBU_PBMC/staroutAligned.sortedByCoord.out.bam | awk 'BEGIN{FS="\t";OFS="\t"} /^[^#]/ { chrom=$1; pos=$2; ref=$4; alts=$5; datas=$10; split(alts,alt, /,/); split(datas, data_fields, /:/); split(data_fields[2], depths, /,/); print chrom, pos, ref, alt[1], depths[1], depths[2]}' > gtcounts.tsv

qsub -V -cwd -b y -j y -N ase -pe smp 2 ../countASE.sh 6_CIRCA_MIBU_PBMC het_subsetSYD-40370185.dedup.realigned.recalibrated.hc.gvcf.gz.bed
qsub -V -cwd -b y -j y -N ase -pe smp 2 ../countASE.sh 7_CIRCA_MIBU_CD4 het_subsetSYD-40370185.dedup.realigned.recalibrated.hc.gvcf.gz.bed
qsub -V -cwd -b y -j y -N ase -pe smp 2 ../countASE.sh 10_COBU_PBMC het_subsetSYD-40341056.dedup.realigned.recalibrated.hc.gvcf.gz.bed
qsub -V -cwd -b y -j y -N ase -pe smp 2 ../countASE.sh 12_JARE_PBMC het_subsetSYD-40341058.dedup.realigned.recalibrated.hc.gvcf.gz.bed
qsub -V -cwd -b y -j y -N ase -pe smp 2 ../countASE.sh 11_ALTR_PBMC het_subsetSYD-40341072.dedup.realigned.recalibrated.hc.gvcf.gz.bed
qsub -V -cwd -b y -j y -N ase -pe smp 2 ../countASE.sh 9_AKAG_CD34 het_subsetSYD-40434090.dedup.realigned.recalibrated.hc.gvcf.gz.bed
qsub -V -cwd -b y -j y -N ase -pe smp 2 ../countASE.sh 8_AKAG_PBMC het_subsetSYD-40434090.dedup.realigned.recalibrated.hc.gvcf.gz.bed
qsub -V -cwd -b y -j y -N ase -pe smp 2 ../countASE.sh GEC_RNA het_subsetGEC.dedup.realigned.recalibrated.hc.gvcf.gz.bed


#some of these bastards took ages for samtools to run so be patient, coz tuncated output doesnt fly



/share/ClusterShare/software/contrib/briglo/bam-readcount_build/bin/bam-readcount -q 13 -l het_subsetSYD-40370185.dedup.realigned.recalibrated.hc.gvcf.gz.bed -f /share/ClusterShare/biodata/contrib/briglo/hg19/hs37d5.fa  briglo/RNA/6_CIRCA_MIBU_PBMC/staroutAligned.sortedByCoord.out.bam


readcounts_het_subset16F00093_counts.bed

samtools mpileup -f /share/ClusterShare/biodata/contrib/briglo/hg19/hs37d5.fa -r 1:14105922-14105922 briglo/RNA/6_CIRCA_MIBU_PBMC/staroutAligned.sortedByCoord.out.bam