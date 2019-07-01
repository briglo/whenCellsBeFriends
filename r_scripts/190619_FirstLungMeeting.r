#Objective Make plots that venessa can use to interpret her data

library(Seurat)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(matlab)
load("r_objects/190619_metadata.rdata")
load("r_objects/190531_integrated.rdata")

#cell type between exp
ggplot(md,aes(x=orig_final.ident,fill=orig.ident)) + geom_bar(position='fill') + geom_hline(yintercept=.5,color='yellow') + ylab("proportion") + coord_flip()

#doublet props
isDoub<-ifelse(md$scrublet_sep>.24 & md$demux_PRB.DBL>.4,"BOTH", ifelse(md$scrublet_sep>.24, "DOUB_scrub", ifelse(md$demux_PRB.DBL>.4,"DOUB_demux",NA)))
##plotting
ggplot(md,aes(x=orig_final.ident,fill=isDoub)) + geom_bar(position='fill') + geom_hline(yintercept=.94,color='yellow') + ylab("proportion") + coord_flip()

#cell type for exp 2
ggplot(md2,aes(x=broad_cell_annotation,fill=orig.ident)) + geom_bar() + geom_hline(yintercept=.5,color='yellow') + ylab("proportion") + coord_flip()

######specific request, expresion patterns of specic genes
GOI<-c("GZMA","GZMB","GZMH","GZMK","GZMM","PFN1","FASL","TNFSF10")
ddat<-t(integrated@assays$SCT@data[rownames(integrated@assays$SCT@data) %in% GOI,])

dddat<-cbind(integrated@meta.data,as(ddat,'matrix'))

ggplot(dddat[!grepl("Cancer",dddat$orig_final.ident),],aes(x=orig.ident,fill=orig.ident,y=GZMA)) + geom_violin() + facet_wrap(~orig_final.ident) + geom_jitter(width=.1) 

###starting on generalizable function
#fetching function
GOI<-c("GZMA","GZMB","GZMH","GZMK","GZMM","PFN1","FASL","TNFSF10")
ddat<-t(integrated@assays$SCT@data[rownames(integrated@assays$SCT@data) %in% GOI,])
keep<-!grepl("Mouse|Doublet|Cancer",integrated@meta.data$orig_final.ident)
##generic plotting function only works on a single gene
#ggplot(integrated@meta.data[keep,],aes(x=orig.ident,fill=orig.ident)) + geom_violin(aes(y=ddat[keep,geneID])) + facet_wrap(~orig_final.ident) + geom_jitter(width=.1,aes(y=ddat[keep,geneID])) + ggtitle(geneID)


#plotting all in one go
adat<-data.frame(id=integrated@meta.data$orig_final.ident,sample=integrated@meta.data$orig.ident,data.frame(as(ddat,'matrix')))
tadta<-melt(adat,id.vars=c('id','sample'))
colnames(tadta)<-c('id','sample','gene','exp')
ggplot(tadta[keep,],aes(y=exp,fill=sample,x=sample)) + geom_violin() + geom_jitter(width=.1,size=.1) + facet_grid(gene~id) 




###proportions of cells with postive expression
tmp<-data.frame(id=integrated@meta.data$orig_final.ident,grp=integrated@meta.data$orig.ident,data.frame(as(ddat>0,'matrix')))

ttmp<-melt(tmp,id.vars=c("id","grp"))
colnames(ttmp)<-c('id','sample','gene','exp_pos')
ggplot(ttmp[keep,],aes(fill=exp_pos,x=sample)) + geom_bar() + facet_grid(gene~id) 
#a more flexible approach 

ttmp[keep,] %>% group_by_all() %>% summarise(n = n()) %>% mutate(prop_cells_expressing=n/sum(n)) %>% filter(exp_pos=="TRUE") %>% complete(sample) %>% ggplot(.,aes(x=gene,fill=sample,y=prop_cells_expressing)) + geom_bar(position='dodge',stat='identity') + facet_wrap(~id)

#as "heatmap" looks like shit
ttmp[keep,] %>% group_by_all() %>% summarise(n = n()) %>% mutate(prop_cells_expressing=n/sum(n)) %>% filter(exp_pos=="TRUE") %>% complete(sample) %>% ggplot(.,aes(x=gene,y=sample,fill=prop_cells_expressing)) + geom_tile() + facet_wrap(~id) + scale_fill_gradient(low='black',high='aquamarine')