############functions#################
#QOL
sys<-function() system('open .')

testInput<-function(rawdata,minGenes=200,minCells=3,totExpr=1e4,GacceptHi=4000,NacceptHi=12000,mitcut=0.15){
    require(Seurat)
    require(Matrix)
    message(paste("mingenes=",minGenes))
    message(paste("mincells=",minCells))
    message(paste("total expr=",totExpr))
    
    print("setting up")
    
    elf5 <- CreateSeuratObject(raw.data=rawdata, min.genes = minGenes, min.cells=minCells,  project = "PyMT", names.field=1) 
    print("annotating mito genes")
    mito.genes <- grep("^mt-", rownames(elf5@data), value = T,ignore.case=T)
    percent.mito <- colSums(elf5@raw.data[mito.genes, ])/colSums(elf5@raw.data)
    elf5 <- AddMetaData(elf5, percent.mito, "percent.mito")
   elf5@meta.data$isElfMouse<-ifelse(grepl("WT",elf5@ident),'WT','Elf5')
    
	  print(table(elf5@ident))
	 message("subsetting data")
    message(paste("gene accept hi=",GacceptHi))
    elf5 <- SubsetData(elf5, subset.name = "nGene", accept.high = GacceptHi)
	 print(table(elf5@ident))
	     message(paste("nUMI accept hi=",NacceptHi))
    elf5 <- SubsetData(elf5, subset.name = "nUMI", accept.high = NacceptHi)
	  print(table(elf5@ident))
	  message(paste("mitoCutoff=",mitcut))
	 elf5<-SubsetData(elf5, subset.name = "percent.mito", accept.high = 0.15)
	  print(table(elf5@ident))
	 elf5<-NormalizeData(elf5,normalization.method = "LogNormalize",scale.factor = 10000)
	 elf5<-ScaleData(object = elf5, vars.to.regress = "nUMI")
	 elf5 <- FindVariableGenes(elf5,  mean.function = ExpMean, dispersion.function = LogVMR,  x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
    print(paste(length(elf5@var.genes),"variable genes"))
    return(elf5)
}

mkPCA<-function(suraObj,pcComp=100){
    message("PCs to compare=",pcComp)
    require(Seurat)
    require(Matrix)
    message("making objects and PCAing, go grab a coffee")
    clus<-suraObj
    clus <- RunPCA(clus, pc.genes=clus@var.genes, pcs.compute = pcComp,do.print=F)
    return(clus)
	PCElbowPlot(clus,num.pc=pcComp)
    }
	 
	 addmetadata<-function(seuratObj,dimvec=1:20) {
	 	#seuratObj@meta.data$isElfMouse<-ifelse(grepl("WT",seuratObj@meta.data$orig.ident),'WT','Elf5')#can be removed *if next time, addded to testinpit
	 	message('cell cycle')
	 	seuratObj<-CellCycleScoring(object = seuratObj, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = FALSE)
	 	message('pca, may or may not complain about convergence')
	 	#seuratObj<<-RunPCA(seuratObj,pcs.compute=100,pc.genes = seuratObj@var.genes, do.print = FALSE) 
	 	message('cluster0')
	 	seuratObj <- FindClusters(seuratObj, dims.use = dimvec, resolution = 0,print.output = FALSE, save.SNN = TRUE)
	 	message('cluster tree')
	 	for (res in seq(0.1,2.8,0.1)) seuratObj <- FindClusters(seuratObj, resolution = res, print.output = FALSE)
	 	message('graphing lies')
	 	seuratObj<- RunTSNE(seuratObj, dims.use = dimvec,do.fast=T)
	 	# message('jackstraw')
# 	 	alldat <- JackStraw(alldat, num.pc = 30, num.replicate = 100, do.print = FALSE)
return(seuratObj)
	 }


 addmetadata_cca<-function(seuratObj,alignmentDim=1:20) {
	 	message('cell cycle')
	 	seuratObj<-CellCycleScoring(object = seuratObj, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = FALSE)
	 	message("aligning subspace in ", max(alignmentDim) ," dimensions. stay cool")
	 	seuratObj <- AlignSubspace(object = seuratObj, reduction.type = "cca", grouping.var = "isElfMouse", 
    dims.align = 1:20)
	message('graphing lies')
	 	seuratObj<- RunTSNE(seuratObj, dims.use = 1:20,reduction.type = "cca.aligned", do.fast=T,check_duplicates = FALSE)
		 message('cluster0')
	 	seuratObj <- FindClusters(seuratObj, dims.use = 1:20,reduction.type = "cca.aligned", resolution = 0,print.output = FALSE, save.SNN = TRUE)
	 	message('cluster tree')
	 	for (res in seq(0.1,2.8,0.1)) seuratObj <- FindClusters(seuratObj, resolution = res, reduction.type = "cca.aligned", print.output = FALSE)
return(seuratObj)
	 }

plotTraj<-function(monoclObj,genevec,filename,metadataDesired,...){
message('depending on your gene list, this might take ages \n you have been warned')
pdf(file=paste0("plots/",filename,".pdf"))
tmp<-monoclObj
tmp <- setOrderingFilter(tmp, genevec)
#print(plot_ordering_genes(tmp))
#ncenters <- length(unlist(sapply(metadataDesired, function(x) unique(as(phenoData(HSMM),'data.frame')[,x]))))
tmp <- reduceDimension(tmp, max_components=2,method = 'DDRTree')#,auto_param_selection = F, nCenter = ncenters)
tmp <- orderCells(tmp, tmp <- orderCells(tmp, reverse=FALSE))
#sapply(c(metadataDesired,"Pseudotime"), function(x) print(plot_cell_trajectory(tmp, color_by=x)))
return(tmp)
dev.off()
}

###########plotting!

compPlot<-function(seuratObj=epidat,ccaObj=epi_cca,var_method="pca",groupVar="isElfMouse") {
p1 <- DimPlot(object = seuratObj, reduction.use = var_method, group.by = groupVar, pt.size = 0.5, 
    do.return = TRUE,plot.title='seurat')
p2 <- DimPlot(object = ccaObj, reduction.use = var_method, group.by = groupVar, pt.size = 0.5, 
    do.return = TRUE,plot.title='cca')
plot_grid(p1, p2)
}	 

plotPie<-function(seuratObj=epidat,primeFactor="isElfMouse",res=0.7){
    x<-table(seuratObj@meta.data[,primeFactor],seuratObj@meta.data[,paste0('res.',as.character(res))])
    x<-x/rowSums(x)
     for(i in 1:ncol(x)) pie(x[,i],main=paste('res=',res,'clus=',colnames(x)[i]))
    plot(1)
return(x)
}

plot_spec_Pie<-function(seuratObj=epidat,primeFactor="isElfMouse",res=0.7,clusters){
    x<-table(seuratObj@meta.data[seuratObj@meta.data[,paste0('res.',as.character(res))] %in% clusters,primeFactor],seuratObj@meta.data[seuratObj@meta.data[,paste0('res.',as.character(res))] %in% clusters , paste0('res.',as.character(res))])
    x<-x/rowSums(x)
    for(i in 1:ncol(x)) pie(x[,i],main=paste('res=',res,'clus=',colnames(x)[i]))
    plot(1)
return(x)
}


	 plotstuff=function(seuratObj){
	 	PCAPlot(seuratObj,group.by='orig.ident')
	 	PCElbowPlot(seuratObj,num.pc=100) 
	 	PCElbowPlot(seuratObj,num.pc=50) 
	 	#JackStrawPlot(alldat, PCs = 1:30)
	 	PCAPlot(seuratObj,group.by='res.0')
	 	PCAPlot(seuratObj,group.by='res.0.1')
	 	PCAPlot(seuratObj,group.by='res.1')
	 	PCAPlot(seuratObj,group.by='res.2.2')
	 	PCAPlot(seuratObj,group.by='isElfMouse')
	 	FeaturePlot(seuratObj,c("percent.mito",'nGene','nUMI','Elf5'),cols.use=c('grey','blue'))
	 	FeaturePlot(seuratObj,c("Epcam",'PyMT','Cd24a','EGFP'),cols.use=c('grey','blue'))
	 	TSNEPlot(seuratObj,group.by='orig.ident')
	 	TSNEPlot(seuratObj,group.by='res.0')
	 	TSNEPlot(seuratObj,group.by='res.0.1')
	 	TSNEPlot(seuratObj,group.by='res.1')
	 	TSNEPlot(seuratObj,group.by='res.2.2')
	 	TSNEPlot(seuratObj,group.by='isElfMouse')
	 }

plo<-function(x="gorsmy_dimRed_x",y="gorsmy_dimRed_y",contVar){
    #print(grep("_x|_y",colnames(epidat@meta.data),value=T))
    #print(grep("HALLMARK",colnames(epidat@meta.data),value=T))
   require(ggplot2)
ggplot(epidat@meta.data,aes_string(x=x , y=y, color=contVar)) + geom_point(alpha=.5) +facet_wrap(~isElfMouse) + scale_colour_gradient2(low = "blue", high = "red") + theme(legend.position="top",legend.text=element_text(size=5))
}

plotMetaScore<-function(seuratObj=epidat,genelist){
    require(Seurat) # for playing with object
    require(ggplot2) # just to make sure can plot
    require(matlab) # for the jet color palette
    require(dplyr) # for the gene contribution transformation
    message(table(genelist %in% rownames(seuratObj@scale.data))['TRUE']," out of ",length(genelist)," genes entered were used to generate score\n") # just shows how many are actually contributing to the score
    hm<-rowSums(seuratObj@scale.data[rownames(seuratObj@scale.data) %in% genelist,]) #
    message("top contributing genes (by percentage) contributing to signature")
    print((hm*100/sum(hm)) %>% sort(decreasing=T) %>% signif(2) %>% head(10) )
    print(ggplot(data.frame(cbind(seuratObj@meta.data[,grep("dimRed",colnames(seuratObj@meta.data))]),metascore=colSums(seuratObj@scale.data[rownames(seuratObj@scale.data) %in% genelist,])),aes(x=gorsmy_dimRed_x,y=gorsmy_dimRed_y,color=metascore)) + geom_point(alpha=.5) + scale_colour_gradientn(colours = jet.colors(7))) #generates the object on the fly and plots it
}


	 
	 ###########enriching
	 getEntrez<-function(seuratobj,geneList='var.genes'){
	     require(biomaRt)
	     require(Seurat)
	 if (geneList=="var.genes") {
	 	return(getBM(attributes=c("ensembl_gene_id","entrezgene","mgi_symbol"),filters="mgi_symbol",values=seuratobj@var.genes,mart=mouse))
	 } else {
	 	return(getBM(attributes=c("ensembl_gene_id","entrezgene","mgi_symbol"),filters="mgi_symbol",values=seuratobj@raw.data,mart=mouse))
	 }
	 }#pulls available data for all variable genesgenes in a seurat object takes a whilew

	 splMar<-function(marker) lapply(split(marker,marker$avg_logFC>0), function(x) mkEnt(split(x$gene,x$cluster))) # splits a FindAllMarkers object into cluster and up/down

	 mkEnt<-function(listOids){
	     require(clusterProfiler)
	     require(ReactomePA)
	     return(lapply(listOids,function(x) as.character(na.omit(hasEntrez$entrezgene[match(x,hasEntrez$mgi_symbol)]))))
	  } # returns a trimmed list of entrez ids from a list of gene symbols

	 reactomeClusts<-function(mkEntObj) {
	     require(clusterProfiler)
	     require(ReactomePA)
	     message("making object, takes a while")
	     return(compareCluster(geneCluster = mkEntObj, fun = "enrichPathway",organism='mouse',universe=as.character(na.omit(hasEntrez$entrezgene)),pvalueCutoff=0.05,readable=T)
	 )
	 } #finds enriched reactome clusters for a list of entrez ids requires a background getEntrez object called hasEntrez

	 hallmarkClusts<-function(mkEntObj) {
	     require(clusterProfiler)
	     require(ReactomePA)
	     message("making object, takes a while")
	     return(compareCluster(geneCluster = mkEntObj, fun = "enricher", organism='mouse', universe=as.character(na.omit(hasEntrez$entrezgene)), pvalueCutoff=0.05, TERM2GENE=hallmark)
	 )
	 }#finds enriched hallmark clusters for a list of entrez ids requires a background getEntrez object called hasEntrez and a hallmark object from hallmark_mouse_TERM2GENE.rdata'


	 keggClusts<-function(mkEntObj) {
	     require(clusterProfiler)
	     require(ReactomePA)
	     message("making object, takes a while")
	     return(compareCluster(geneCluster = mkEntObj, fun = "enrichKEGG",organism='mouse',universe=as.character(na.omit(hasEntrez$entrezgene)),pvalueCutoff=0.05)
)
}
