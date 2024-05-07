### SCRIPT FOR IMMUNE GENE SIGNATURES (E.G. CANCER IMMUNITY CYCLE GENE SIGNATURES) ###

##Set Directories and Load FSA Package##

.libPaths(c(""))# Path to libraries
.libPaths() 
setwd("")# Set working directory
knitr::opts_chunk$set(echo = TRUE)

library(FSA)
set.seed(123)

##Load gene signatures##
gene_sig <- read.delim2("gene_signatures_entrezid.txt", header=FALSE, row.names=1) # gene signatures in entrezID
# Add any set of gene signatures in entrezID format

##Load Counts Data and remove negative values##

# Transcript counts data post-normalisation and sequencing batch correction using ComBat
df <- read.csv(".csv", row.names=1) # can also import in .txt format using read.table

df[df< 0] <- 0 #Remove negative values (script cannot work with these)

##Set up functions##

inject.dots <- function(df) {names(df) <- sub(" ", ".", names(df));df}

resetZeros<-function(expression_data, to){
  newdata<-expression_data
  newdata[newdata==0] <- to
  return (newdata)
}
create_empty_table <- function(num_rows, num_cols) {
  frame <- data.frame(matrix(NA, nrow = num_rows, ncol = num_cols))
  return(frame)
}


getMinVal<-function(expression_data){
  min<-0
  for (i in seq(1,length(expression_data[1,]))) {
    for (cur in expression_data[,i]) {
      if (min == 0){ min=cur}
      if (cur>0) {
        if (cur<min){
          min<-cur
        }
      }
    }	
  }
  return (min)
}

check_sig<-function(signatures,expression_data){
  
  samples<-colnames(expression_data)
  sigs<-row.names(signatures)
  res <- create_empty_table(length(sigs),4)
  colnames(res)<-c("size","percentageCovered","missingInstances", "existingInstances")
  rownames(res)<-sigs
  sample<-samples[1]
  
  for (i in sigs){
    siggenes<-signatures[i,]
    sigvec<-unlist(strsplit(as.vector(siggenes[1]),","))
    vec2<-c()
    vec3<-c()
    for (si in sigvec){ #si means each gene
      l<-expression_data[si,sample] #give me for each sample of the matrix the gene
      if (!is.finite(l)){
        vec3<-c(vec3,si)
      } else {
        vec2<-c(vec2,si)
      }
      
    }
    
    amount_there<-length(vec2)
    res[i,1]<-length(sigvec)
    res[i,2]<-amount_there/length(sigvec)
    o<-paste(vec3, sep=",", collapse=",")
    o2<-paste(vec2, sep=",", collapse=",")
    res[i,3]<-o
    res[i,4]<-o2
    
  }
  return (inject.dots(res))
}

get_sig_gm<-function(signatures,expression_data){
  
  samples<-colnames(expression_data)
  sigs<-row.names(signatures)
  res <- create_empty_table(length(samples),length(sigs))
  rownames(res)<-samples
  colnames(res)<-sigs
  for (sample in samples){
    for (i in sigs){
      siggenes<-signatures[i,]
      sigvec<-unlist(strsplit(as.vector(siggenes[1]),","))
      vec2<-c()
      for (si in sigvec){
        l<-expression_data[si,sample]
        vec2<-c(vec2,l)
      }
      gm<-geomean(vec2, na.rm=TRUE) #change to mean if needed? 
      res[sample,i]<-gm
      
    }
  }
  return (inject.dots(res))
}

get_sig_m<-function(signatures,expression_data){
  
  samples<-colnames(expression_data)
  sigs<-row.names(signatures)
  res <- create_empty_table(length(samples),length(sigs))
  rownames(res)<-samples
  colnames(res)<-sigs
  for (sample in samples){
    for (i in sigs){
      siggenes<-signatures[i,]
      sigvec<-unlist(strsplit(as.vector(siggenes[1]),","))
      vec2<-c()
      for (si in sigvec){
        l<-expression_data[si,sample]
        vec2<-c(vec2,l)
      }
      gm<-mean(vec2, na.rm=TRUE) #change to mean if needed? 
      res[sample,i]<-gm
      
    }
  }
  return (inject.dots(res))
}

##Set up expression data##

expression_data<-as.data.frame(df) #should be a df, not a matrix !
min_val<-getMinVal(expression_data) #does min values have an effect? 
min_val<-min_val/10
expression_data<-resetZeros(expression_data, min_val)

mysigs_gene_sig<-gene_sig
m<-check_sig(mysigs_gene_sig,expression_data)

##Calculate geometric mean using ENTREZID genes and save output##
n_gene_sig_gm<-get_sig_gm(mysigs_gene_sig,expression_data)
write.csv(n_gene_sig_gm,"XXX.csv") # write.table for saving as .txt
