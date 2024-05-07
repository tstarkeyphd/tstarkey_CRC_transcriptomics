### RNAseq COUNT CRC CLASSIFICATION (CMS AND CRIS) ###

# Set libpaths
.libPaths(c(""))

# Set working directory
setwd("")
knitr::opts_chunk$set(echo = TRUE)

##Consensus molecular subtype classification##
library(CMSclassifier)
set.seed(123)

data<-read.csv("FILENAME.csv", header=T, row.names=1)# RNAseq counts post-normalisation and batch adjustment (if applicable)
Rfcms <- CMSclassifier::classifyCMS(data,method="RF")[[3]]
write.csv(Rfcms,"FILENAME.csv")
rm(data,Rfcms)

# Figure generation
library(pheatmap)
library(heatmap3)

data<-read.csv(".csv", header=T, row.names=1)

breaksList = seq(-350, 300, by = 50)

cc1<-palette(c(
  "steelblue4","steelblue3","steelblue2","steelblue1","white",
  "goldenrod1","goldenrod2","goldenrod3","goldenrod4"))
palette()
                           
cc2<-palette(c(
  "palevioletred4","palevioletred3","palevioletred2","palevioletred1","white",
  "goldenrod1","goldenrod2","goldenrod3","goldenrod4"))
palette()
                                                          
cc3<-palette(c(
  "seagreen4","seagreen3","seagreen2","seagreen1","white",
  "goldenrod1","goldenrod2","goldenrod3","goldenrod4"))
palette()

cc4<-palette(c(
  "palevioletred4","palevioletred3","palevioletred2","palevioletred1",
  "white","steelblue1","steelblue2","steelblue3","steelblue4"))
palette()

cc5<-palette(c(
  "seagreen4","seagreen3","seagreen2","seagreen1","white",
  "steelblue1","steelblue2","steelblue3","steelblue4"))
palette()

cc6<-palette(c(
  "seagreen4","seagreen3","seagreen2","seagreen1","white",
  "palevioletred1","palevioletred2","palevioletred3","palevioletred4"))
palette()

pheatmap(
  data, cluster_rows = FALSE, cluster_cols = FALSE,
  clustering_method = "complete",clustering_distance_cols = "euclidean",# Only important if cluster_rows/cols=TRUE
  color = colorRampPalette(cc1)(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same length of breaksList)
  breaks = breaksList) # Sets the breaks of the color scale as in breaksList
# Change cc colour palette (cc1) as per requirements

##Colorectal intrinsic subtype classification##
library(CRISclassifier)

data<-read.csv(".csv", header=T)# RNAseq counts post-normalisation and batch adjustment (if applicable)
cris_classifier(input.exp.filename = data, output.name="cris", nresmpl=1)
rm(data)