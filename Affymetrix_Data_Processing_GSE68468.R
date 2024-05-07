### AFFYMETRIX PROCESSING CaArray_notte-00422 ###


# Set libpaths
.libPaths(c(""))

# Set working directory
setwd("")
knitr::opts_chunk$set(echo = TRUE)

# Load packages
library(limma)
library(oligo)
library(Biobase)
library(affy)
library(GEOquery)
library(gcrma)

# Load data
downloadedAffyFiles <- list.files(path = "", pattern = "CEL.gz$",full.names=TRUE)#Path to files

AffyData <- ReadAffy(filenames = downloadedAffyFiles)

## Background correction and normalisatio (2 methods using gcrma)

set.seed(123)

# Option 1 - rma() function (used for my analyses)
rma<-rma(AffyData, background=TRUE, normalize=TRUE, target="core") # can change core to probeset
write.csv(rma,"AffyData_rma_core.csv")

# Option 2 - gcrma() function
gcrma <- gcrma(AffyData, background=TRUE, normalize=TRUE, target="core") # can change core to probeset
write.csv(gcrma,"AffyData_gcrma_core.csv")


# Create table with affymetrix probe ids and corresponding gene symboles
require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(
  mart=mart,
  attributes=c(
    "affy_primeview",
    "ensembl_gene_id",
    "gene_biotype",
    "external_gene_name"),
  filter = "affy_primeview",
  values = rownames(exprs(rma)), uniqueRows=TRUE) # can change filter and rma to 'gcrma'

write.csv(annotLookup,"FILE.csv")

# Average values for genes with multiple probe ids (limma)
data <- read.csv("AffyData_rma_core.csv", header=T)

data[,1]-> gene #take SYMBOL genes
as.character(gene) -> gene

data_summarised <- limma::avereps(data, ID = gene)

data2 <-as.data.frame(data_summarised)

write.csv(data2,"AffyData_Gene_Counts.csv")
