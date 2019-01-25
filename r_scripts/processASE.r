
library(GenomicRanges)
library(parallel)
setwd('/share/ScratchGeneral/circa/briglo/ASE')
fnames<-dir(pattern="tsv")
snames<-gsub("_refAltCounts.tsv","",fnames)

rdat<-mclapply(fnames, function(X) read.table(X, header=F, stringsAsFactors=F),mc.cores=length(fnames))
names(rdat)<-snames
lapply(rdat,dim)
trdat<-mclapply(rdat, function(x) return(x[x$V5+x$V6>20,]),mc.cores=length(fnames))
lapply(trdat,dim)
gtrdat<-mclapply(trdat, function(x) return (GRanges(seqnames=x$V1,ranges=IRanges(start=x$V2-1,end=x$V2),mcols=data.frame(ref=x$V3,alt=x$V4,refcount=x$V5, altcount=x$V6))),mc.cores=length(fnames))

save(gtrdat,file="ASE_readCountsGRanges.rdata")

#binomialTestFromMark
pfunk = function (ad, dp) min(1, 2*pbinom(min(ad, dp-ad), dp, 0.5, lower.tail=TRUE))
pv<-lapply(gtrdat, function(x) apply(cbind(x$mcols.refcount,x$mcols.refcount+x$mcols.altcount),1,function(y) pfunk(y[1],y[2])))
pa<-lapply(pv, function(x) p.adjust(x,method="BH"))
for (i in 1:length(pa)) gtrdat[[i]]$pdj_pval=pa[[i]]
save(gtrdat, file='ASE_readCountsGRanges_binomTested.rdata')

#GetExons affected by ase
tase<-lapply(gtrdat, function(x) x[x$pdj_pval<0.001])
load("~/share/ScratchGeneral/circa/annotations/genes/hg19_ens.ano.rdata")
exons<-hg19_ens[hg19_ens$type=="exon"]

exonsaffectedbyase<-lapply(tase, function(x) exons[exons %over% x])
save(exonsaffectedbyase, file="ASE.rdata")

############
library(GenomicRanges)
library(parallel)
setwd('/share/ScratchGeneral/circa/briglo/ASE')
load("ASE_readCountsGRanges.rdata")
load("/share/ScratchGeneral/circa/annotations/genes/hg19_ens.ano.rdata")

imp<-read.table("briglo/ASE/ImprintedGenes.tsv",header=T,skip=1,stringsAsFactors=F,sep='\t')
impSplit<-split(exons,exons$gene_name %in% imp$GeneID)


spread<-lapply(gtrdat, function(x) {
    tmp<-x$mcols.refcount/rowSums(cbind(x$mcols.refcount,x$mcols.altcount))
    isimp<- x %over% impSplit[[2]]
    mcols(x)$refProp=tmp
    mcols(x)$isImprinted=isimp
    return(x)

})

trimGencode<-function(ids) sub("^(.*)[.].*", "\\1",ids)

tmp<-do.call(rbind,lapply(spread,as.data.frame))
tmp$sample=trimGencode(rownames(tmp))
tmp<-tmp[!tmp$mcols.alt=="<*>",]
tmp$pbio<-apply(data.frame(v1=tmp$mcols.refcount,v2=tmp$mcols.refcount+tmp$mcols.altcount), 1, function(x) pfunk(x[1],x[2])
tmp$adj_pbio<-p.adjust(tmp$pbio,method="BH")


ggplot(tmp[tmp$seqnames %in% c(seq(1,23,1),"X"),],aes(x=refProp,color=seqnames,fill=seqnames)) + geom_density(alpha=0.2) + facet_wrap(~sample)
ggsave("Spreads_byChr_bySamp.pdf") #looks like a cluster on chr6

ggplot(tmp[tmp$seqnames %in% c(seq(1,23,1),"X"),],aes(x=refProp,color=isImprinted,fill=isImprinted)) + geom_density(alpha=0.2) + facet_wrap(~sample)
ggsave("Spreads_byImprinted_bySamp.pdf")

ggplot(tmp[tmp$seqnames %in% c(seq(1,23,1),"X"),],aes(x=refProp,color=sample,fill=sample)) + geom_density(alpha=0.2) + facet_wrap(~seqnames)
ggsave("Spreads_bySamp_byChrom.pdf")


ggplot(tmp[tmp$seqnames %in% c("X") & tmp$adj_pbio<0.05,],aes(x=start,y=abs(0.5-refProp),color=sample)) + geom_point() + facet_wrap(~sample)
ggsave("Spreads_bychrom_bySamp_MANHATTAN.pdf")


pfunk = function (ad, dp) min(1, 2*pbinom(min(ad, dp-ad), dp, 0.5, lower.tail=TRUE))


tmp$pbio<-apply(data.frame(v1=tmp$mcols.refcount,v2=tmp$mcols.refcount+tmp$mcols.altcount), 1, function(x) pfunk(x[1],x[2])
tmp$adj_pbio<-p.adjust(tmp$pbio,method="BH")


#makig into a GRanges again
ASE<-GRanges(seqnames=tmp$seqnames, ranges=IRanges(start=tmp$start,end=tmp$end))
colnames(tmp)<-gsub("mcols.","",colnames(tmp))
mcols(ASE)<-tmp[,c(6:14)]


load("annotations/genes/hg19_ens.ano.rdata")
exons<-hg19_ens[hg19_ens$type=='exon']

maxASE<-ASE[ASE$adj_pbio<0.001]
smaxASE<-split(maxASE,maxASE$sample)
smaxASEGenes<-lapply(smaxASE, function(x) table(exons$gene_name[exons %over% x]))
lapply(smaxASEGenes, function(x) head(x[order(x)]))
#shows sense, lacks rigour, try another way


#Annotating ASE with which exon(s) it hits
x<-findOverlaps(ASE,exons)
geneNames<-exons$gene_name[subjectHits(x)]
sgeneNames<-split(geneNames,queryHits(x))
test<-unlist(lapply(sgeneNames, function(x) paste(unique(x),collapse=";")))
dftest<-as.data.frame(test,stringsAsFactors=F)
names(ASE)<-1:length(ASE)
ASE$geneExonHit<-dftest[names(ASE),1]


#make a summary of ASE candidates vs the actual numner of ASE cases per exon target. This may be faulty
oi<-table(ASE$geneExonHit[ASE$adj_pbio<0.001],ASE$sample[ASE$adj_pbio<0.001])
oi2<-table(ASE$geneExonHit,ASE$sample)
m<-match(rownames(oi),rownames(oi2))
toi2<-oi2[m,]

library(reshape2)
moi<-melt(oi)
mtoi2<-melt(toi2)
chu<-data.frame(cbind(moi,mtoi2))
tchu=chu[!chu$value==0 & !chu$value.1==0,]
ASE_byGeneExonHit<-tchu[,c(1,2,3,6)]
colnames(ASE_byGeneExonHit)<-c('gene','sample','numberASEadjp0.001','totalCandForASE')
ASE_byGeneExonHit$prop<-ASE_byGeneExonHit$numberASEadjp0.001/ASE_byGeneExonHit$totalCandForASE

plotCands<-function(totCutoff=5,propcutoff=.5,geneomit="none"){
ggplot(ASE_byGeneExonHit[ASE_byGeneExonHit$totalCandForASE>totCutoff & ASE_byGeneExonHit$prop>propcutoff & !grepl(geneomit,ASE_byGeneExonHit$gene),],aes(x=gene,y=prop,fill=sample)) + geom_bar(stat='identity') + facet_wrap(~sample,scales='free') + coord_flip()
}

save(tmp,ASE,ASE_byGeneExonHit,file='ASEtableForPlotting.rdata')
