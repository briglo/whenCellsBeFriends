# whenCellsBeFriends
in silico identification and deconvolution of scRNA doublets

## comments
* This is a mix of things, dont freak out, there should be a master process below
* Rationalised R functions in /R Recently updated for Seurat3,
* related to github.com/briglo/scFuncs
* R objects that were particularly painful to get are in  /share/ClusterShare/thingamajigs/SCCG/briglo/venchi_doub

## a simple (if time comsuming) Seurat workflow
```R
setwd("PATH/TO/DATA")
datadirs<- paste0(dir(),"/outs/filtered_feature_bc_matrix")
bysample<-makeSeuratList(datadirs)
integrated<-buildMasterSeurat(bysample)

# look at population composition 
plotSankey(integrated,c("orig.ident","Phase"))

# look at reactome enrichment for clusters
markers<-makeReactomePipe(integrated,"SCT_snn_res.0.8")
for (i in 1:length(markers$CP_result)) dotplot(markers$CP_result[[i]]) + ggtitle(names(markers$CP_result)[i])

# look at combined expression of some marker genes in UMAP space
integrated<-makeMetaScore(integrated,markers$markers$gene[1:50])
ggplot(integrated@meta.data,aes(x=UMAP_1,y=UMAP_2,color=metascore)) + geom_point(alpha=.5) + scale_colour_gradientn(colours = jet.colors(7))
```

# questions
1) Are there more intelligent ways to identify doublets
  * probably not by me- all the buzzwords have been used
  * however semisupervised using model based on similarity scores to things like scmca maybe

2) Can we distinguish in situ doublets from accidental doublets
  * based on cellphone interactivity score maybe
3) Can we use in-situ doublets to artifically reconstruct a spatial map
  * Probably some awful hidden markoff model based on 2

  # plan
  1) Try to get allele aware count data (start with Kang et al- should be an intermediate step of demuxlet)
  2) Run demuxlet and scrublet over current data
  3) Compare each output to known inputs and intersect

![justification](https://github.com/briglo/whenCellsBeFriends/blob/master/images/190509_justification.png)

# speculation
  1) Intersection demuxlet and scrublet should yield accidental duplicates (demux only), falsly "interacting" duplicates (both) and possible interactors(scrublet only)
  2) These may give insight as to the determinants of cell- interactions that are retained through fluidics (and my guess is that solid tissue based libraries should have a higher proportion of interactors)
3) Wonder if cluster based might be the way to go? (define clusters of interest, use teichman to ID dups between them, extract bam, run freebayes ahh i might have hit my glucose intake for the day) for future application in an unknown space... would shoehorn better in to existing pipelines.


# background
1) Doublets occur when two cells are encased in the microfluidics of a scRNA device
2) A crude method of identifying doublets is by number of transcripts 
3) Doublets have weird transcriptomes and are typically discarded from analysis

## relevant reading
### Kang Nat Biotech 2018 Demuxlet
"harnesses genetic variation to determine the genetic identity of each droplet... and identify driples containing cells from differenet individuals".  
  *made a module, runs fine*


###  Zheng Nature Comms 2016
- Poisson Dist estimate of multiplet rate ~ 9%
- cite that dissecting similar celltype multiplet will be hard

**method to get genotype aware counts (could not find specific method) prefers freeBayes SNV detection over more standard approaches (3% MAF)**

### scrublet wolock biorxiv-python
- KNN estimate of doublet from simulating doublets from count matrix, better than nGene)
- (plus nice commentary on most approaches to detect multiplets)
- cannot detect a pair if the "parent" population is underepresented (cites example of macrophage & erythroid cell)
- calculates phiD which can estimate the potential impact of "doublets" in the data

### doubletDecon Depasquale biorxiv- R
like scrublet but it takes into account identified cell states (by seurat blah blah, trims, reclusters, rescues possible transitional states...)
No support on github since August last year, no documnetation on file formats, not example data.

### doubleFinder biorxiv 2018 
like scrublet but a poorer written paper

### Stoeckius Nat Methods 2017
Spikein species at low copy: ~4% doublet

### Bloom Peer J (?) 2018
- Estimation of multiplet frequency using cell mixing experiements (will accept >2)
- Clumping would change this

### Ilicic 2016 Genome Biol
ML (SVM) for low qual, can ID multiples but is dependent on "closeness" to training set...  Not currently interested




