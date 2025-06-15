##################################################################################################################
### Description: Calculate tumor mes/adrn scores from bulk RNA sequencing data
###              

library(dplyr)
library('gtable')
library(data.table)
library(parallel)
library(tidyr)
workdir = '/rocker-build/gmkf_nbl_somatic/'

#import clinical data
clinData = read.csv(paste0(workdir,'Data/tumor_clinical_data.tsv'),sep='\t',header=TRUE,row.names=NULL,stringsAsFactors=FALSE)

#import RNA readcount data
expMat = read.table(paste0(workdir,'Data/RNAseq_rsem_coding_TPM.tsv'),sep='\t',stringsAsFactors=FALSE,header=TRUE)
gNames = rownames(expMat)
sNames = colnames(expMat)
expMat = data.table::transpose(expMat)
colnames(expMat) = gNames
rownames(expMat) = sNames

#list of coding genes
codGenes = readRDS(paste0(workdir,'References/hg38_coding_genes.RDS'))

#filter expression matrix by coding genes only
expMatCod = expMat[,colnames(expMat) %in% codGenes]

##ADRN vs MES signature
#import adrn and mes genes
amGenes = read.table(paste0(workdir,'References/adrn_mes_genes.csv'),sep=',',stringsAsFactors=FALSE)
mGenes = amGenes$V1[amGenes$V2 == 'MES']
aGenes = amGenes$V1[amGenes$V2 == 'ADRN']
amGenes = setNames(list(mGenes,aGenes),nm=c('MES','ADRN'))

###calculate percentile rank for each gene in each sample
##for each sample, calculate percentile rank for each gene

allGenes = colnames(expMatCod)
mGenes = mGenes[mGenes %in% colnames(expMatCod)]
aGenes = aGenes[aGenes %in% colnames(expMatCod)]
totGenes = length(allGenes)

#convert genes to index number in column order
mIX = as.numeric(sapply(mGenes,FUN=function(x) { c(1:totGenes)[allGenes == x] }))
aIX = as.numeric(sapply(aGenes,FUN=function(x) { c(1:totGenes)[allGenes == x] }))

#calculate mesenchymal and adrn scores
maScores = sapply(c(1:nrow(expMatCod)),FUN=function(x) {
  vec = as.numeric(expMatCod[x,])
  vecOrd = order(vec)
  
  mRanks = as.numeric(mclapply(as.list(mIX),FUN = function(y){
    placeInVec = c(1:totGenes)[vecOrd == y]
    rank = placeInVec/totGenes
    return(rank)
  },mc.cores = 8))
  
  aRanks = as.numeric(mclapply(as.list(aIX),FUN = function(z){
    placeInVec = c(1:totGenes)[vecOrd == z]
    rank = placeInVec/totGenes
    return(rank)
  },mc.cores = 8))
  
  mAvg = mean(mRanks)
  aAvg = mean(aRanks)
  return(c(mAvg,aAvg))
})

colnames(maScores) = rownames(expMatCod)
rownames(maScores) = c('MES_SCORE','ADRN_SCORE')
maScores = t(maScores)

maDF = as.data.frame(maScores)

#add A/M ratio
maDF$AM_RATIO = maDF$ADRN_SCORE/maDF$MES_SCORE

#add color for plotting
maDF$CLASS = sapply(maDF$AM_RATIO,FUN=function(x) {
  if (x<=0.9) {
    return('Mesenchymal')
  } else if (x>0.9 & x <=1) {
    return('Slight Mesenchymal')
  } else if (x>1 & x <=1.1) {
    return('Slight Adrenergic')
  } else {
    return('Adrenergic')
  }
})

maDF$COLOR = sapply(maDF$CLASS,FUN=function(x) {
  if (x=='Mesenchymal') {
    return('Green')
  } else if (x=='Slight Mesenchymal') {
    return('Gray')
  } else if (x=='Slight Adrenergic') {
    return('Gray')
  } else {
    return('Red')
  }
})

#plot
par(mar=c(5,5,5,5))
plot(maDF[,2],maDF[,1],xlim=c(0.5,1),ylim=c(0.5,1),bg=maDF$COLOR,xlab='ADRN_SCORE',ylab='MES_SCORE',cex=3.5,cex.lab=2,cex.axis=2,pch=21)
abline(a=c(0,1))
abline(a=c(0.08,1),lty=2)
abline(a=c(-0.069,1),lty=2)

#confirm values of MES and ADRN scores
total_adrn = sum(maDF$AM_RATIO > 1)