BKKD DEseq2
================

``` r
GeneCountData <- read.delim("all.combine.cntTable", header = TRUE, sep = "\t", row.names = 1)
GeneCountDataInt <- as.matrix(GeneCountData)
storage.mode(GeneCountDataInt) = "integer"
GeneCountDataInt2<-GeneCountDataInt[,c(1,3,4,6)]
Design <- read.delim("comsample.txt", header = TRUE, sep = "\t", row.names = 1)
```

Run Deseq2
==========

``` r
library("DESeq2")

##generating a DESeq dataset matrix object using the counts and design input files
gdds <- DESeqDataSetFromMatrix(countData = GeneCountDataInt2, colData = Design, design = ~ condition)

##performing DE analysis
gdds_de <- DESeq(gdds)

##generating normalized counts for each gene or transcript for plots
dds_rld <- rlog(gdds)

##generating a results object from DE analysis (see above)
dds_res <- results(gdds_de)

##generating subset of DE results for individual comparisons (for me, the "genotype" design, comparing MVKKO vs MVWT; could also be something like "tissue" design comparing liver vs kidney or liver vs testis [whatever is labeled in your file]")
res <- results(gdds_de, contrast=c("condition","A","B"),alpha=0.05)
```

MA plot \#plotMA(res, ylim=c(-5,5))
===================================

Label the gene
==============

``` r
wantgene <- row.names(res)[which(res$padj < 0.05 & res$log2FoldChange > 1)]
wantgene
```

    ## [1] "ENSDARG00000009922" "ENSDARG00000028148" "ENSDARG00000031222"
    ## [4] "ENSDARG00000035610" "ENSDARG00000074253" "ENSDARG00000096563"

``` r
plotMA(res, ylim=c(-5,5))
with(res[wantgene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, wantgene, pos=2, col="dodgerblue")
})
```

![](Deseq2_files/figure-markdown_github/unnamed-chunk-3-1.png)
