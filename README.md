# whenCellsBeFriends
in silico identification and deconvolution of scRNA doublets

# background
1) Doublets occur when two cells are encased in the microfluidics of a scRNA device
2) A crude method of identifying doublets is by number of transcripts 
3) Doublets have weird transcriptomes and are typically discarded from analysis


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


## comments
* This is a mix of things, dont freak out, there should be a master process below
* uses functions from github.com/briglo/scFuncs
* R objects that were particularly painful to get are in  /share/ClusterShare/thingamajigs/SCCG/briglo/venchi_doub

##  where i got to
### 1) Trimming bam files for "useful" records
```bash
for i in * ; do qsub -pe smp 4 -b y -j y -N rmd -V -cwd "samtools view -h "$i"/outs/possorted_genome_bam.bam hg19_1 hg19_2 hg19_3 hg19_4 hg19_5 hg19_6 hg19_7 hg19_8 hg19_9 hg19_10 hg19_11 hg19_12 hg19_13 hg19_14 hg19_15 hg19_16 hg19_17 hg19_18 hg19_19 hg19_20 hg19_21 hg19_22 | samtools view -h -S -q 10 -F 3844 - | samtools sort -@ 4 -m 4G -o "$i"/outs/human_trim_sort_"$i".bam -" ; done
for i in * ; do qsub -b y -j y -N rmd -V -cwd samtools rmdup -s --reference /share/ScratchGeneral/briglo/scRNA/chrcha/anno/renamed_hg19.fa "$i"/outs/human_trim_sort_"$i".bam "$i"/outs/rmdup_human_trim_sort_"$i".bam ; done 
for i in * ; do qsub b y -j y -N rmd -V -cwd samtools index "$i"/outs/rmdup_human_trim_sort_"$i".bam ; done 
```

### 2) extracting "pure" barcode lists for bam subsetting
```R
md<-read.csv("extraData/pilot_1_metadata.csv",header=T,stringsAsFactors=F)
x<-split(md,md$orig.ident)
extBC<-function(metaData){
bcs<- gsub("_[0-9]",'-1',metaData$Barcode)
return(list('cancer'=bcs[grep("Cancer",metaData$final.ident)],
        'normal'=bcs[metaData$final.ident %in% c("Monocyte","CD14+,Macrophage","CD4,T-cell,2","T-Reg","CD8,T-cell,2","Dendritic,cells","NK","CD16+,Macrophage","pDC","CD8,T-cell,1","B-cell","CD4,T-cell,1","Plasmablast")]
        ))
}
bclist<-lapply(x,extBC)
write.table(bclist$Sample_A$cancer,quote=F,col.names=F,row.names=F,file="data/hg19/VENCHI_SampleACITE/outs/cancer_barcodes.txt")
write.table(bclist$Sample_A$normal,quote=F,col.names=F,row.names=F,file="data/hg19/VENCHI_SampleACITE/outs/normal_barcodes.txt")
write.table(bclist$Sample_B$cancer,quote=F,col.names=F,row.names=F,file="data/hg19/VENCHI_SampleBCITE/outs/cancer_barcodes.txt")
write.table(bclist$Sample_B$normal,quote=F,col.names=F,row.names=F,file="data/hg19/VENCHI_SampleBCITE/outs/normal_barcodes.txt")
```

### 3) extracting bam records
```bash
module load briglo/subset-bam/1.0_precompiled
cd /share/ScratchGeneral/briglo/scRNA/venchi/data/hg19/VENCHI_SampleACITE/outs
qsub -V -cwd -pe smp 4 -b y -j y subset-bam --bam rmdup_human_trim_sort_VENCHI_SampleACITE.bam --bam-tag CB:Z --cell-barcodes normal_barcodes.txt --cores 4 --out-bam rmdup_human_trim_sort_VENCHI_SampleACITE_normal.bam
qsub -V -cwd -pe smp 4 -b y -j y subset-bam --bam rmdup_human_trim_sort_VENCHI_SampleACITE.bam --bam-tag CB:Z --cell-barcodes cancer_barcodes.txt --cores 4 --out-bam rmdup_human_trim_sort_VENCHI_SampleACITE_cancer.bam
samtools index rmdup_human_trim_sort_VENCHI_SampleACITE_normal.bam
samtools index rmdup_human_trim_sort_VENCHI_SampleACITE_cancer.bam
cd /share/ScratchGeneral/briglo/scRNA/venchi/data/hg19/VENCHI_SampleBCITE/outs
qsub -V -cwd -pe smp 4 -b y -j y subset-bam --bam rmdup_human_trim_sort_VENCHI_SampleBCITE.bam --bam-tag CB:Z --cell-barcodes normal_barcodes.txt --cores 4 --out-bam rmdup_human_trim_sort_VENCHI_SampleBCITE_normal.bam
qsub -V -cwd -pe smp 4 -b y -j y subset-bam --bam rmdup_human_trim_sort_VENCHI_SampleBCITE.bam --bam-tag CB:Z --cell-barcodes cancer_barcodes.txt --cores 4 --out-bam rmdup_human_trim_sort_VENCHI_SampleBCITE_cancer.bam
samtools index rmdup_human_trim_sort_VENCHI_SampleBCITE_normal.bam
samtools index rmdup_human_trim_sort_VENCHI_SampleBCITE_cancer.bam
```

### 4) Executing FreeBayes and running Demuxlet: Two approaches, I personally had a harder time interpreting the output from the ReadGroup approach
#### The method i prefer (per file of cancer/normal, from sc_split)
```bash
 samtools view -S -b -q 10 -F 3844 original.bam > target.bam
 # run freebayes
 module load briglo/freebayes/v1.2.0-4-gd15209e 
  freebayes -f /share/ScratchGeneral/briglo/scRNA/chrcha/anno/renamed_hg19.fa -iXu -C 2 -q 1 target.bam > snv.vcf
 #QC trim
  module load briglo/bcftools/1.9
bcftools filter -i 'QUAL>30' -O z -o snv.vcf snv_filter.vcf

#combine vcfs into a single VCF for demux
echo "SampleBCITE_cancer" > can
 bcftools reheader --samples can -o SampleBCITE_cancer_filter_rename.vcf.gz SampleBCITE_cancer_filter.vcf.gz
 tabix -p vcf SampleBCITE_normal_filter.vcf.gz
  tabix -p vcf SampleBCITE_cancer_filter_reheader.vcf.gz
  bcftools merge -O z -o SampleBCITE_merged.vcf.gz SampleBCITE_normal_filter.vcf.gz SampleBCITE_cancer_filter_rename.vcf.gz

#run demuxlet  
module load briglo/demuxlet/1.9
 demuxlet --sam rmdup_human_trim_sort_VENCHI_SampleACITE.bam  --vcf snv_filter.vcf --out ./outputs/demuxlet/Sample  --field GT --geno-error 0.001
```

#### the readgroup method (freebayes_2.sh)
```bash
# merge cancer/normal bam files after renanaming readgroups
module load briglo/freebayes/v1.2.0-4-gd15209e 
cd /share/ScratchGeneral/briglo/scRNA/venchi/data/hg19/VENCHI_SampleACITE
 qsub -V -cwd -pe smp 2 -N renMerge_A -m e -M b.gloss@garvan.org.au -b y -j y "bamaddrg -b rmdup_human_trim_sort_VENCHI_SampleACITE_cancer.bam -s ACITE_cancer -b rmdup_human_trim_sort_VENCHI_SampleACITE_normal.bam -s ACITE_normal | samtools sort -o ACITE_MergeRG.bam -"
   cd ../../VENCHI_SampleBCITE/outs
   qsub -V -cwd -pe smp 2 -N renMerge_B -m e -M b.gloss@garvan.org.au -b y -j y "bamaddrg -b rmdup_human_trim_sort_VENCHI_SampleBCITE_cancer.bam -s BCITE_cancer -b rmdup_human_trim_sort_VENCHI_SampleBCITE_normal.bam -s BCITE_normal | samtools sort -o BCITE_MergeRG.bam -"

#running freebayes scatter-gather coz i am impatient
hdir="/share/ScratchGeneral/briglo/scRNA/venchi"
module load briglo/freebayes/v1.2.0-4-gd15209e 
for i in  `ls data/hg19/` ; do mkdir "$hdir"/outputs/freebayes/"$i"/ 
for j in `cat "$hdir"/source/chrom.txt` ; do \
qsub -pe smp 2 -b y -j y -N fb_"$i"_"$j" -m e -M b.gloss@garvan.org.au -V -cwd "samtools view -h "$hdir"/data/hg19/"$i"/outs/*CITE_MergeRG.bam "$j" | freebayes -f /share/ScratchGeneral/briglo/scRNA/chrcha/anno/renamed_hg19.fa -iXu -C 2 -q 1 --min-coverage 12 --stdin  -v "$hdir"/outputs/freebayes/"$i"/"$i"_"$j"_normal_SNV.vcf"
done
done

###compress
mv *gz byChr_old/
mv *tbi byChr_old/
module load briglo/samtools/1.5
for i in *vcf ; do qsub -V -cwd -b y -j y bgzip -i "$i" ; done
 rm bgzip.o*
###combine
module load pethum/vcftools/gcc-4.4.6/0.1.15
ls *normal*.gz > norm.txt
vcf-concat -f norm.txt > ACITE.vcf.gz
mkdir byChr
mv VENCHI_SampleBCITE_hg19_* byChr
mv norm.txt byChr

#QC
module load briglo/bcftools/1.9
bcftools filter -i 'QUAL>30' -O z -o ACITE_filter.vcf.gz ACITE.vcf.gz
tabix -p vcf ACITE_filter.vcf.gz

 # run demuxlet, had to use GT as opposed to GP (does imputation calc this? yes but probably not nec, important to give alpha if possible(seyhan), playing with alpha gives a hot mess)
  module load briglo/demuxlet/1.9
  
  qsub -V -cwd -pe smp 2 -b y -j y -m e -M b.gloss@garvan.org.au -N dm demuxlet --sam ./data/hg19/VENCHI_SampleACITE/outs/rmdup_human_trim_sort_VENCHI_SampleACITE.bam --vcf ./outputs/freebayes/VENCHI_SampleACITE/ACITE_filter.vcf.gz --out ./outputs/demuxlet/SampleACITE_RG  --field GT --geno-error 0.001 
  
  #--alpha 0.5 --alpha 0.3
  
  qsub -V -cwd -pe smp 2 -b y -j y -m e -M b.gloss@garvan.org.au -N dm demuxlet --sam ./data/hg19/VENCHI_SampleBCITE/outs/rmdup_human_trim_sort_VENCHI_SampleBCITE.bam --vcf ./outputs/freebayes/VENCHI_SampleACITE/ACITE_filter.vcf.gz --out ./outputs/demuxlet/SampleBCITE_RG  --field GT --geno-error 0.001 
  
  #--alpha 0.5 --alpha 0.4
```

### 5) compare in R
```R
fnam<-dir(pattern=PREFIXOFCHOICE)
  dm<-lapply(fnam, function(x) data.frame(read.table(x,header=T,stringsAsFactors=F),row.names=1)[-1,]) 
  md<-read.csv("../../extraData/pilot_1_metadata.csv",header=T,stringsAsFactors=F)
  rownames(md)<-md$Barcode

lapply(dm, function(X) table(X$ALPHA))
smoothScatter(dm[[1]]$PRB.DBL,dm[[1]]$PRB.SNG1)
smoothScatter(dm[[2]]$PRB.DBL,dm[[2]]$PRB.SNG1)

rownames(dm[[1]])<-gsub("-[0-9]","_1",rownames(dm[[1]]))
rownames(dm[[2]])<-gsub("-[0-9]","_2",rownames(dm[[2]]))
cdm<-data.frame(do.call(rbind,dm))
cdm<-cdm[rownames(md),]
plot(cdm$PRB.DBL,cdm$PRB.SNG1)

#combining var (cdm_var) alpha (0.3 and 0.6 based on ~reads/transcript ratio between cancer and epithelials) and usual (0.5, cdm)
cdm$isDoub<-ifelse(cdm$PRB.DBL>0.6 & cdm$PRB.SNG1<.95,"doublet",'singlet')
cdm_var$isDoub<-ifelse(cdm_var$PRB.DBL>0.55 & cdm_var$PRB.SNG1<.95 & grepl("_2",colnames(cdm_var)),"doublet",ifelse(cdm_var$PRB.DBL>0.6 & cdm_var$PRB.SNG1<.95,'doublet','singlet')
colnames(cdm_var)<-paste0("varAlpha_",colnames(cdm_var))
allcdm<-data.frame(cbind(cdm,cdm_var))
save(allcdm,file="../../r_objects/190606_demux_alpha_RG_comparison.rdata")
```



# history
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
