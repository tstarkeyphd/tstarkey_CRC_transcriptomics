### RNAseq COUNT PROCESSING, NORMALISATION AND BATCH EFFECT CORRECTION ###

# Set libpaths
.libPaths(c(""))

# Set working directory
setwd("")
knitr::opts_chunk$set(echo = TRUE)



# Change gene SYMBOL to ENTREZID

library(org.Hs.eg.db)
library(Homo.sapiens)
library(AnnotationDbi)
library(GenomicFeatures)
library(stringi)
library(digest)
library(XML)
library(curl)
library(devtools)
library(Biobase)

# Load data
df <- read.csv(".csv", header=T) #HTSeq-processed RNAseq counts data

# Transform gene symbols to entrezID (required for downstream analyses)
df[,1]-> geneid #take SYMBOL genes
as.character(geneid) -> geneid

genes_ENTREZ<-mapIds(
  Homo.sapiens,keys=geneid, column="ENTREZID", 
  keytype ="SYMBOL") #Change to Entrezid

cbind(df, genes_ENTREZ) -> df2
df2 <- na.omit(df2) #remove nas - a lot bc of symbol + ".number"
df2<-df2[!duplicated(df2$genes_ENTREZ),] #none ? 
rownames(df2) <- df2$genes_ENTREZ
dim(df2)
df2 <- df2[,-142]
dim(df2)
df2->df

df<-df[,-(464),drop=FALSE] 

write.csv(df,".csv")
rm(df2,genes_ENTREZ)


#Data normalisation and batch correction using ComBat
library(car)
library(edgeR)
library(survival)
library(sva)

set.seed(123)

#log2 +1 transformation for normalisation 
df<-read.csv(".csv", header=T, row.names = 1)
df[df< 0] <- 0
df<-as.matrix(t(df))

df_norm<- cpm(df, log = TRUE, prior.count = 1) 
df_norm<-as.data.frame(t(df_norm))

#remove batch effect using ComBat
batch<-read.csv(".csv",header=T)
cbind(df_norm,batch)->df_norm
nsamples <- ncol(df_norm)

colnames(df_norm) -> batch
as.data.frame(batch) -> batch
rownames(batch) <- colnames(df_norm)
as.character(batch$batch) -> batch
matr <- as.matrix(df_norm)

combat <- ComBat(
  dat =  matr, batch=batch, prior.plots=TRUE,
  par.prior=TRUE) #uses the parametric Bayesian adjustments
combat<-as.data.frame(combat)

write.csv(combat,".csv")