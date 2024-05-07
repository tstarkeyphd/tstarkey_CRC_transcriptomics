### DESeq2 Differential Expression Analysis Script ###


# Set directories
.libPaths(c(""))
.libPaths() 
setwd("")
knitr::opts_chunk$set(echo = TRUE)



#Differential expression analysis using DESeq2


#Load DESeq2
library(DESeq2)

getwd() # Find the working directory from file pane and past into line below
directory <- "Path/To/Working/Directory"

# directory should be the location where the sampleFiles (.out) are located

# can merge individual sample files (i.e. ctrl1.counts, ctrl2.counts, etc.)
sampleFiles <- grep(".out", list.files(directory), value=T)#RNAseq counts from HTSeq

# view sampleFiles
sampleFiles

# can designate different batches of samples (i.e. different sequencers,
# PE vs SE, different library preps (eg. BATCH1 vs BATCH2))

# set sampleConditions and sampleTable for experimental conditions within the excel file sampleTable.csv

# view sampleTable
sampleTable <- read.csv("sampleTable.csv")
head(sampleTable)



###Pathway 1 - DESeq2 by Wald Test###


# Run first step of DESeq2
ddsHTseq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design = ~condition) 

##  Set condition levels
colData(ddsHTseq)$condition <- factor(colData(ddsHTseq)$condition,
                                      levels = c('A','B'))
# gut of DESeq2 analysis
dds <-DESeq(ddsHTseq, test="Wald")
res <- results(dds)
# order results by padj value (most significant to least)
res <- res[order(res$padj), ]
head(res)

# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj 
# save data results and normalized reads to csv!
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized=T)), by='row.names',sort=F)
names(resdata)[1] <- 'gene'
head(resdata)
write.csv(resdata, file="FILENAME_Wald.csv")



###Pathway 2 - DESeq2 by Likelihood Ratio Test (LRT)###
#Useful for accounting for sequencing batches or assessing a second condition

# Run first step of DESeq2
ddsHTseq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design = ~condition + Batch) 

##  Set condition levels
colData(ddsHTseq)$condition <- factor(colData(ddsHTseq)$condition,
                                      levels = c('A','B'))
# gut of DESeq2 analysis
dds <-DESeq(ddsHTseq, test="LRT", reduced=~Batch)
res <- results(dds)
# order results by padj value (most significant to least)
res <- res[order(res$padj), ]
head(res)

# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj 
# save data results and normalized reads to csv!
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized=T)), by='row.names',sort=F)
names(resdata)[1] <- 'gene'
head(resdata)
write.csv(resdata, file="FILENAME_LRT.csv")


### VOLCANO PLOTS ###

res <- read.csv("volcano.csv", header=TRUE)#DESeq2 data (differential gene expression)
head(res)

# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Differential Expression ", xlim=c(-0.6,0.6)),frame=FALSE)

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>0.25), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, pvalue<.05 & abs(log2FoldChange)>0.5), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
# p-value and log2 fold change cut-offs can be modified

# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(res, pvalue<.05 & abs(log2FoldChange)>5), textxy(log2FoldChange, -log10(pvalue), labs=gene, cex=.7))
# p-value and log2 fold change cut-offs can be modified
