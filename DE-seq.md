DE-seq on 8 samples
================

``` r
setwd("/Volumes/Junybug/BKKD_Junybug/2nd/")
GeneCountData <- read.delim("allgene.fc.count", sep = "\t", row.names = 1)
```

*processes file to make sure that it is in the right format*

``` r
GeneCountDataInt <- as.matrix(GeneCountData)
storage.mode(GeneCountDataInt) = "integer"
```

GeneCountDataInt2&lt;-GeneCountDataInt\[,c(1,3,4,6)\]
\#GeneCountDataInt3&lt;-GeneCountDataInt\[,c(1,4,5,6,8,9)\]

*Loading in your design file*

``` r
Design <- read.delim("sample.txt", header = TRUE, sep = "\t", row.names = 1)
```

*Loading “DEseq2” software*

``` r
library("DESeq2")
```

*generating a DESeq dataset matrix object using the counts and design
input files*

``` r
gdds <- DESeqDataSetFromMatrix(countData = GeneCountDataInt, colData = Design, design = ~ condition + time)
```

*performing DE analysis*

``` r
gdds_de <- DESeq(gdds)
```

*generating normalized counts for each gene or transcript for plots*

``` r
dds_rld <- rlog(gdds)
```

*plots PCA of samples using the groups specified in your design file
(for me, genotype and serum)*

``` r
plotPCA(dds_rld, intgroup=c("condition", "time"))
```

res &lt;- results(gdds\_de,
contrast=c(“condition”,“Exp”,“Sc”),alpha=0.05) \#\#tells you total \# of
genes DE at adj. p &lt; 0.1 level summary(res) plotMA(res, ylim=c(-5,5))

write.csv(res, file=“testwolalal.csv”)

**It looks good. It separates conditions and time.** <br> *now lets do
8hpf only*

**8hpf**

``` r
GeneCountDataInt2<-GeneCountDataInt[,c(1,4,5,8)]
head(GeneCountDataInt2)
```

``` r
Design2 <-Design[c(1,4,5,8),]
```

``` r
gdds2 <- DESeqDataSetFromMatrix(countData = GeneCountDataInt2, colData = Design2, design = ~ condition)

gdds_de2 <- DESeq(gdds2)

dds_rld2 <- rlog(gdds2)

plotPCA(dds_rld2, intgroup=c("condition"))
```

*generating a results object from DE analysis*

``` r
#dds_res2 <- results(gdds_de2)
res2 <- results(gdds_de2, contrast=c("condition","Exp","Sc"),alpha=0.05)
##tells you total # of genes DE at adj. p < 0.1 level
summary(res2)
```

*MA plt*

``` r
plotMA(res2, ylim=c(-5,5))
```

\#\#writing DE results to CSV file write.csv(res2,
file=“KD\_2nd\_fc\_0126.csv”)

**12hpf**

``` r
GeneCountDataInt3<-GeneCountDataInt[,c(2,3,6,7)]
head(GeneCountDataInt3)
```

``` r
Design3 <-Design[c(2,3,6,7),]
```

``` r
gdds3 <- DESeqDataSetFromMatrix(countData = GeneCountDataInt3, colData = Design3, design = ~ condition)

gdds_de3 <- DESeq(gdds3)

dds_rld3 <- rlog(gdds3)

plotPCA(dds_rld3, intgroup=c("condition"))
```

*generating a results object from DE analysis*

``` r
#dds_res2 <- results(gdds_de2)
res3 <- results(gdds_de3, contrast=c("condition","Exp","Sc"),alpha=0.05)
##tells you total # of genes DE at adj. p < 0.1 level
summary(res3)
```

*MA plt*

``` r
plotMA(res3, ylim=c(-5,5))
```

\#\#writing DE results to CSV file write.csv(res3,
file=“KD\_2nd\_12hpf\_0126.csv”)

*combine 2 RNA-seq*

``` r
GeneCountData0<- read.delim("../all.combine.cntTable", header = TRUE, sep = "\t", row.names = 1)
GeneCountDataInt0 <- as.matrix(GeneCountData0)
storage.mode(GeneCountDataInt0) = "integer"

head(GeneCountDataInt0)
head(GeneCountDataInt2)
```

``` r
#library(tidyverse)

GeneCountDataInttest <- GeneCountDataInt0 %>% merge(GeneCountDataInt2, by=0, all=TRUE)
head(GeneCountDataInttest)
write.csv(GeneCountDataInttest, file="combine.csv")

test <- read.delim("combine.csv", header = TRUE, sep = ",", row.names = 1)
test <- as.matrix(test)
storage.mode(test) = "integer"
```

``` r
Designtest <- read.delim("sample_combine.txt", header = TRUE, sep = "\t", row.names = 1)
```

``` r
gddstest <- DESeqDataSetFromMatrix(countData = test, colData = Designtest, design = ~ condition)

gdds_detest <- DESeq(gddstest)

dds_rldtest <- rlog(gddstest)

plotPCA(dds_rldtest, intgroup=c("condition"))
```

``` r
test1 <-test[,c(1,3,4,6,7,8,9,10)]
```

``` r
Designtest1 <- read.delim("sample_combine2.txt", header = TRUE, sep = "\t", row.names = 1)
```

``` r
gddstest1 <- DESeqDataSetFromMatrix(countData = test1, colData = Designtest1, design = ~ condition)

gdds_detest1 <- DESeq(gddstest1)

dds_rldtest1 <- rlog(gddstest1)

plotPCA(dds_rldtest1, intgroup=c("condition"))
```

``` r
#dds_res2 <- results(gdds_de2)
restest1 <- results(gdds_detest1, contrast=c("condition","Exp","Sc"),alpha=0.05)
##tells you total # of genes DE at adj. p < 0.1 level
summary(restest1)
```

``` r
plotMA(restest1, ylim=c(-5,5))
```

write.csv(restest1, file=“testtest\_0126.csv”)

**heatmap**

``` r
library("pheatmap")
sampleDists <- dist(t(assay(dds_rld2)))
library("RColorBrewer")
sampleDists <- dist(t(assay(dds_rld2)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(dds_rld2$replicate, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```

*my heatmap* library(“pheatmap”) sampleDists &lt;-
dist(t(assay(dds\_rld))) library(“RColorBrewer”) sampleDists &lt;-
dist(t(assay(dds\_rld))) sampleDistMatrix &lt;- as.matrix(sampleDists)
rownames(sampleDistMatrix) &lt;-
paste(dds\_rld*r**e**p**l**i**c**a**t**e*, *v**s**d*type, sep=“-”)
colnames(sampleDistMatrix) &lt;- NULL colors &lt;- colorRampPalette(
rev(brewer.pal(9, “Blues”)) )(255) pheatmap(sampleDistMatrix,
clustering\_distance\_rows=sampleDists,
clustering\_distance\_cols=sampleDists, col=colors)

\*MA plot \#plotMA(res, ylim=c(-5,5)) wantgene &lt;-
row.names(res)\[which(res$padj &lt; 0.05 & res$log2FoldChange &gt; 1)\]
wantgene wantgene2 &lt;-
row.names(res)\[which(res$padj &lt; 0.05 & res$log2FoldChange &lt; -1)\]
wantgene2 wantgene3 &lt;- c(wantgene, wantgene2)

wantgene\_manual &lt;- c(“ENSDARG00000028148”, “ENSDARG00000031222”,
“ENSDARG00000074253”, “BHIKHARI\_I-int:LTR:LTR”,
“BHIKHARI\_LTR:LTR:LTR”) names(wantgene\_manual)

plotMA(res, ylim=c(-4,4)) with(res\[wantgene\_manual, \], {
points(baseMean, log2FoldChange, col=“grey”, cex=2, lwd=2)

})

text(baseMean, log2FoldChange, text(baseMean, log2FoldChange,
wantgene\_manual, pos=2, col=“dodgerblue”), pos=2, col=“dodgerblue”)

\*loading in necessary software to generate sample distance heatmaps
library(“RColorBrewer”) library(“gplots”)

\#\#generating objects/variables necessary for sample distance heatmaps
\#\#note that most of these do not need to be changed, with the
exception of your appropriate group specified in your design file \[for
me, “genotype”\]) distsStageRL &lt;- dist(t(assay(dds\_rld)))
gene\_stage\_mat &lt;- as.matrix(distsStageRL)
rownames(gene\_stage\_mat) &lt;- colnames(gene\_stage\_mat) &lt;-
with(colData(gdds),paste(genotype, sep=" : “)) g\_stage\_col &lt;-
colorRampPalette(brewer.pal(9,”GnBu"))(100) hc1 &lt;-
hclust(distsStageRL)

*generating heatmap* heatmap.2(gene\_stage\_mat,
Rowv=as.dendrogram(hc1), symm=TRUE, trace=“none”, col =
rev(g\_stage\_col), margin=c(13,13))

*generating subset of DE results for individual comparisons (for me, the
“genotype” design, comparing MVKKO vs MVWT; could also be something like
“tissue” design comparing liver vs kidney or liver vs testis \[whatever
is labeled in your file\]")* res &lt;- results(gdds\_de,
contrast=c(“condition”,“A”,“B”))
