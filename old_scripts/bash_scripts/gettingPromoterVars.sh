####in R
library(rtracklayer)
load("~/share/ScratchGeneral/circa/annotations/genes/hg19_ens.ano.rdata")
tx<-hg19_ens[hg19_ens$type=="transcript"]
nam<-tx$transcript_id
tss<-ifelse(as.character(strand(tx))=="+",start(tx),end(tx))
upstream<-ifelse(as.character(strand(tx))=="+",tss-400,tss+400)
downstream<-ifelse(as.character(strand(tx))=="+",tss+250,tss-250)

starts<-apply(cbind(tss,upstream), 1, min)
ends<-apply(cbind(tss,upstream), 1, max)
tx_prom_upstream<-GRanges(seqnames=seqnames(tx), ranges=IRanges(start=starts,end=ends),strand=as.character(strand(tx)),names=nam,gene_name=tx$gene_name)

starts<-apply(cbind(downstream,upstream), 1, min)
ends<-apply(cbind(downstream,upstream), 1, max)
tx_prom<-GRanges(seqnames=seqnames(tx), ranges=IRanges(start=starts,end=ends),strand=as.character(strand(tx)),names=nam,gene_name=tx$gene_name)

save(tx_prom,tx_prom_upstream,file="~/share/ScratchGeneral/circa/briglo/promvars/ens_hg19_promoterregions.rdata")

export(reduce(tx_prom),file="/share/ScratchGeneral/circa/annotations/genes/hg19_ens_transcript_promoters_reduce.bed","bed")

#########################

module load pethum/bedtools/gcc-4.4.6/2.25.0
for i in `cat tmp` ; do \
qsub -V -cwd -b y -j y -N btp_"$i" -pe smp 1 \
"bedtools intersect -wa -header \
-a vcfs/"$i" \
-b annotations/genes/hg19_ens_transcript_promoters_reduce.bed > ./promSubset"$i"" 

for i in promSubset* ; do cat "$i" | awk 'BEGIN {FS="\t";OFS="\t"} ($7 == "PASS" && ($10  ~ /^0\/1:/)) {print $1,$2-1,$2,$4,$5}' > het_promSubset"$i" ; done