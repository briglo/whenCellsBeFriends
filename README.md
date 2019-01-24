# whenCellsBeFriends
in silico identification and deconvolution of scRNA doublets

# background
1) Doublets occur when two cells are encased in the microfluidics of a scRNA device
2) A crude method of identifying doublets is by number of transcripts 
3) Doublets have weird transcriptomes and are typically discarded from analysis

## relevant reading
### Kang Nat Biotech 2018 Demuxlet
"harnesses genetic variation to determine the geneitc identity of each droplet... and identify driples containing cells from differenet individuals".  *made a module, runs fine*


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

### doubleFinder biorxiv 2018 
like scrublet but a poorer written paper

### Stoeckius Nat Methods 2017
Spikein species at low copy: ~4% doublet

### Bloom Peer J (?) 2018
- Estimation of multiplet frequency using cell mixing experiements (will accept >2)
- Clumping would change this

### Ilicic 2016 Genome Biol
ML (SVM) for low qual, can ID multiples but is dependent on "closeness" to training set...  Not currently interested




# questions
1) Are there more intelligent ways to identify doublets
  * probably not by me- all the buzzwords have been used
  * however semisupervised using model based on similarity scores to things like scmca maybe

2) Can we distinguish in situ doublets from accidental doublets
  * based on cellphone interactivity score maybe
3) Can we use in-situ doublets to artifically reconstruct a spatial map
  * Probably some awful hidden markoff model based on 2

