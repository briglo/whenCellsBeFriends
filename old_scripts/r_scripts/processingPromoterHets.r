library(GenomicRanges)
library(parallel)
setwd('/share/ScratchGeneral/circa/briglo/promvars')
####generating GRanges of results from gettingPromoterVars.sh####
fnames<-dir(pattern="het_promSubsetpromSubset")
snames<-gsub("het_promSubsetpromSubset|.hc.vqsr.vcf.gz","",fnames)
rdat<-mclapply(fnames, function(X) read.table(X, header=F, stringsAsFactors=F),mc.cores=length(fnames))
 lapply(rdat,dim)
gtrdat<-mclapply(rdat, function(x) return (GRanges(seqnames=x$V1,ranges=IRanges(start=x$V2,end=x$V3),mcols=data.frame(ref=x$V4,alt=x$V5))),mc.cores=length(fnames))
names(gtrdat)<-c('JARE','ALTR','COBU','MIBU','AKAG') # based on sample spreadsheet
save(gtrdat,file="promoter_hets_GRanges.rdata")
#########




load("~/share/ScratchGeneral/circa/briglo/promvars/promoter_hets_GRanges.rdata")
load("~/share/ScratchGeneral/circa/briglo/promvars/ens_hg19_promoterregions.rdata")

promaffectedbyhet<-lapply(gtrdat, function(x) tx_prom[tx_prom %over% x])
purepromaffectedbyhet<-lapply(gtrdat, function(x) tx_prom_upstream[tx_prom_upstream %over% x])
save(promaffectedbyhet,purepromaffectedbyhet,file="transcriptPromotersAffectedbyhet_GRanges.rdata")


load("../ASE/transcriptExonsAffectedbyASE.rdata")


####e.g. ranks by numbers
table(purepromaffectedbyhet$'JARE'$names %in% exonsaffectedbyase$'12_JARE_PBMC'$transcript_id) #15% seems high probably coz "promoter" region is in the exons and thus counted... anyway as a handle on what's what did following
 genes<-table(purepromaffectedbyhet$'JARE'$gene_name[purepromaffectedbyhet$'JARE'$names %in% exonsaffectedbyase$'12_JARE_PBMC'$transcript_id])
 genes<-genes[order(genes,decreasing=T)]
 head(genes,20)

 #### for all genes... not ranked
 gene_ids_ASE_assoc_with_prom_hets=list(
     "COBU_PBMC"=unique(exonsaffectedbyase$'10_COBU_PBMC'$gene_id[exonsaffectedbyase$'10_COBU_PBMC'$transcript_id %in% purepromaffectedbyhet$'COBU'$names]),
     "ALTR_PBMC"=unique(exonsaffectedbyase$'11_ALTR_PBMC'$gene_id[exonsaffectedbyase$'11_ALTR_PBMC'$transcript_id %in% purepromaffectedbyhet$'ALTR'$names]),
     "JARE_PBMC"=unique(exonsaffectedbyase$'12_JARE_PBMC'$gene_id[exonsaffectedbyase$'12_JARE_PBMC'$transcript_id %in% purepromaffectedbyhet$'JARE'$names]),
     "MIBU_PBMC"=unique(exonsaffectedbyase$'6_CIRCA_MIBU_PBMC'$gene_id[exonsaffectedbyase$'6_CIRCA_MIBU_PBMC'$transcript_id %in% purepromaffectedbyhet$'MIBU'$names]),
     "MIBU_CD4"=unique(exonsaffectedbyase$'7_CIRCA_MIBU_CD4'$gene_id[exonsaffectedbyase$'7_CIRCA_MIBU_CD4'$transcript_id %in% purepromaffectedbyhet$'MIBU'$names]),
     "AKAG_PBMC"=unique(exonsaffectedbyase$'8_AKAG_PBMC'$gene_id[exonsaffectedbyase$'8_AKAG_PBMC'$transcript_id %in% purepromaffectedbyhet$'AKAG'$names]),
     "AKAG_CD34"=unique(exonsaffectedbyase$'9_AKAG_CD34'$gene_id[exonsaffectedbyase$'9_AKAG_CD34'$transcript_id %in% purepromaffectedbyhet$'AKAG'$names])
 )

#   PDE4DIP     SEPT2    GIGYF2      IL32     CARD8     FCHO1     WDR41     NAA60
#        52        35        29        26        24        23        23        21
#    NAP1L1 LINC00969     CPNE1    NPEPPS   POLR2J3     RBM39      RNF4   NADSYN1
#        21        20        18        18        18        18        18        17
#     ZNF83  CALCOCO2      MDM2    MRPL52
#        17        16        16        16



###to do, exclude transcribed variants from crossover