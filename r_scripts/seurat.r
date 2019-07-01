screen
qrsh -pe smp 12
module load briglo/R/3.6.0
source .profile
cd /share/ScratchGeneral/briglo/scRNA/venchi
R

library(Seurat)
library(Matrix)
library(reticulate)
use_python("/share/ClusterShare/software/contrib/briglo/miniconda3/envs/magic/bin/python")
source("source/seuratFunctions.r")
load("r_objects/cell.cyclegenes.rdata")
snam<-dir("data/hg19")

sink("maybeouts.txt")
for(i in 1:length(snam)){

id<- lapply(snam, function(x) {
    ddir<-paste0(getwd(),"/data/hg19/",x,"/outs/filtered_feature_bc_matrix/")
#system(paste0("gunzip ",ddir,"*"))
rd <- Read10X(ddir)
rd<-rd[grepl("hg19",rownames(rd)),]
rownames(rd)<-gsub("hg19_","",rownames(rd))
return(CreateSeuratObject(counts=rd,project=gsub('VENCHI_Sample',"",x), min.cells=3,min.features=200))
})


anchors <- FindIntegrationAnchors(object.list = id, dims = 1:30)
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
integrated@active.assay="RNA"
integrated[["percent.mt"]] <- PercentageFeatureSet(integrated, pattern = "^MT-")
integrated<-subset(integrated, subset = percent.mt < 10)
integrated<-CellCycleScoring(object = integrated, s.features = toupper(s.genes), g2m.features = toupper(g2m.genes), set.ident = FALSE)
integrated@active.assay="integrated"
integrated<-SCTransform(integrated, vars.to.regress = c("percent.mt","nFeature_RNA","S.Score","G2M.Score"), verbose = FALSE)


integrated<-RunPCA(integrated, verbose = FALSE)
ElbowPlot(integrated)
integrated<-RunUMAP(integrated, dims = 1:30, verbose = T)

DimPlot(integrated, reduction = "umap", pt.size = 0.1) ; dev.off()

integrated<-FindNeighbors(integrated, dims = 1:30, verbose = FALSE)
integrated<-FindClusters(integrated, verbose = FALSE)


integrated<-RunTSNE(integrated, genes.use=integrated@var.genes, dims.use=1:17, do.fast=T)

md<-data.frame(cbind(integrated@reductions$umap@cell.embeddings[,1:2], integrated@reductions$pca@cell.embeddings[,1:3], integrated@reductions$tsne@cell.embeddings[,1:2],integrated@meta.data))
z<-read.csv("extraData/pilot_1_metadata.csv",header=T,stringsAsFactors=F)
rownames(z)<-z$Barcode # i ignored mouse transcripts (might want to not? for accidental buddy rate dep on when added)
tz<-z[rownames(md),]
colnames(tz)<-paste0("orig_",colnames(tz))
md<-data.frame(cbind(md,tz[,c(8,23)]))

integrated@meta.data<-md
save(integrated,file='r_objects/190516_CCAintegrateAllCells.rdata')

##### scrublet
writeMM(integrated@assays$RNA@counts,file="outputs/seurat/combined.human.mtx")
write.table(rownames(integrated@assays$RNA@counts),file="outputs/seurat/genes.tsv",quote=F,row.names=F,col.names=F)
write.table(colnames(integrated@assays$RNA@counts),file="outputs/seurat/barcodes.tsv",quote=F,row.names=F,col.names=F)
#scrublet.py
 scr<-read.table('outputs/scrublet/Sample_BCITE_doubletScores_seuratCombined_HUMAN.tsv')
 A<-read.table('outputs/scrublet/Sample_ACITE_doubletScores_filteredMat.tsv')
ab<-read.table("data/hg19/VENCHI_SampleACITE/outs/filtered_feature_bc_matrix/barcodes.tsv")
rownames(A)<-gsub("-[1-9]","_1",ab$V1)

B<-read.table('outputs/scrublet/Sample_BCITE_doubletScores_filteredMat_ALL.tsv')
bb<-read.table("data/hg19/VENCHI_SampleBCITE/outs/filtered_feature_bc_matrix/barcodes.tsv")
rownames(B)<-gsub("-[1-9]","_2",bb$V1)

allb<-read.table("outputs/seurat/barcodes.tsv")
rownames(scr)<-allb$V1

sdm<-data.frame(rbind(A,B))[rownames(scr),]
colnames(scr)<-"scrublet_comb"
scr$scrublet_sep<-sdm

#### demuxlet
#freebayes.sh
dmA<-data.frame(read.table("outputs/demuxlet/SampleACITE.best",header=T)[-1,],row.names=1,stringsAsFactors=F)
dmB<-data.frame(read.table("outputs/demuxlet/SampleBCITE.best",header=T)[-1,],row.names=1,stringsAsFactors=F)
rownames(dmA)<-gsub("-[1-9]","_1",rownames(dmA))
rownames(dmB)<-gsub("-[1-9]","_2",rownames(dmB))
cdm<-data.frame(rbind(dmA,dmB))[rownames(integrated@meta.data),]
colnames(cdm)<-paste0("demux_",colnames(cdm))

#using readgroup renaming pipeline in freebayes_2 (ALPHA SHOULD BE 0.5 as of JUNE coz fusking with alpha was wrong)
# fdmA_RG<-data.frame(read.table("outputs/demuxlet/SampleACITE_RG.best",header=T)[-1,],row.names=1,stringsAsFactors=F)
# dmB_RG<-data.frame(read.table("outputs/demuxlet/SampleBCITE_RG.best",header=T)[-1,],row.names=1,stringsAsFactors=F)
# rownames(dmA_RG)<-gsub("-[1-9]","_1",rownames(dmA_RG))
# rownames(dmB_RG)<-gsub("-[1-9]","_2",rownames(dmB_RG))
# cdm_RG<-data.frame(rbind(dmA_RG,dmB_RG))[rownames(integrated@meta.data),]
# colnames(cdm_RG)<-paste0("demux_RG_",colnames(cdm_RG))
### RG didnt really do much at all except confuse things with multiallelic stuff. For another time
##results live in 190606_demux_alpha_RG_comparison.rdata

plot(cdm$demux_PRB.DBL,cdm_RG$demux_RG_PRB.DBL)
integrated@meta.data<-data.frame(cbind(integrated@meta.data,cdm_RG[,c(4,11,12,20,21)])
save(integrated, file="r_objects/190531_integrated.rdata")

