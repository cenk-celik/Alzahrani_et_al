# Introduction

This repository describes step-by-step RNA sequencing (RNA-seq) analysis conducted in human PRE1 cell line by Alzahrani _et al._ (in progress), where ER stress sensor IRE1 is knocked out (IRE1_KO) and compared with untreated RPE1 cells (IRE1_WT).

# RNA-seq steps
## Step 1: Obtaining count matrices

We employed [```salmon```](https://salmon.readthedocs.io/en/latest/) (v 1.9.0) for obtaining the count matrices from the ```_.fastq_``` files using [Gencode Release 42](https://www.gencodegenes.org/human/) with parameters ```--gcBias``` and ```--validateMappings```. We also generated ```decoys.txt``` as described in ```salmon``` [vignette](https://salmon.readthedocs.io/en/latest/salmon.html?highlight=decoy#preparing-transcriptome-indices-mapping-based-mode). The ```salmon``` script for obtaining count matrices can be found [here](/scripts/quantifier.sh).

## Step 2: ```tximport``` for ```DESeq2```

After quantification, we used [```tximport```](https://bioconductor.org/packages/release/bioc/html/tximport.html) (v 1.26.1) to import counts generated at step 1. For genome-wide annotation, we used [Gencode Annotation Release 42](https://www.gencodegenes.org/human/). See the workflow [here](scripts/salmon2tximport.R).

## Step 3: ```DESeq2``` object and global response

Using the ```txi``` object that includes all count matrices, we generated ```dds``` object using [```DESeq2``` package](http://master.bioconductor.org/packages/release/bioc/html/DESeq2.html) (v 1.38.3). Following filtering out genes with low counts, we visualised PCA plot, top 500 variable genes and sample distances. For the script, click [here](scripts/DESeq2overall.r).

## Step 4: Differentially expressed genes

Using this [script](scripts/DESeq2results.R), we observed differentially expressed genes (DEGs) in IRE1_KO, compared with control.

## Step 5: Gene Ontology analysis

For over-representation and gene set enrichment analysis (GSEA), we used [```clusterProfiler```](http://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) (v 4.6.0), [```enrichplot```](http://bioconductor.org/packages/release/bioc/html/enrichplot.html) (v 1.18.3) and [```DOSE```](http://bioconductor.org/packages/release/bioc/html/DOSE.html) (v 3.24.2) packages for [Gene Ontology](scripts/clusterProfiler.R) analyses, as well as [```ReactomePA```](http://bioconductor.org/packages/release/bioc/html/ReactomePA.html) (v 1.42.0) for [Reactome](scripts/reactome.R).

## Step 6: KEGG pathway analysis

Similar to Gene Ontology analysis, we also looked at [```KEGG```](https://www.kegg.jp) pathways. The script can be found [here](scripts/kegg.R).
