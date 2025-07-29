##################################################################################################################
### Description: This script finds structural variant mutations that potentially disrupt coding genes
###              - Helper functions are imported from the script sv_overlap_functions.R that perform intersection
###                with gene bed files
###              

library(parallel)
library(data.table)
library(dplyr)
workdir = '/rocker-build/gmkf_nbl_somatic/'
source(paste0(workdir,'Scripts/sv_overlap_functions.R'))


###import bed file of hg38 transcripts
bed = read.table(paste0(workdir,'References/hg38_ucsc_transcripts.bed'),sep='\t',
                 stringsAsFactors=FALSE,header=TRUE,row.names=NULL)
bedEx = read.table(paste0(workdir,'References/hg38_ucsc_exons.bed'),sep='\t',
                   stringsAsFactors=FALSE,header=TRUE,row.names=NULL)
chrs = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
         'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
         'chr18','chr19','chr20','chr21','chr22','chrX','chrY')
bed$Chromosome = factor(bed$Chromosome,levels=chrs)
bedEx$Chromosome = factor(bedEx$Chromosome,levels=chrs)
#filter bed file for coding genes
genes = readRDS(paste0(workdir,'References/hg38_coding_genes.RDS'))
bed = bed %>% filter(Gene %in% genes)
#split bed by chromosome
bed = split(bed,f=bed$Chromosome)
#split exons df into list by gene name
bedGene = split(bedEx,f=bedEx$Gene)


###load example SV filtered events data
#Lumpy
lFiles = list.files(paste0(workdir,'Raw Data/Lumpy Filtered Events Examples'))
lFiles = as.list(sapply(lFiles,FUN=function(x){strsplit(x,split='_')[[1]][1]},USE.NAMES=FALSE))
lsv = mclapply(lFiles,FUN=function(x){
  fName = paste0(workdir,'Raw Data/Lumpy Filtered Events Examples/',x,'_filtered_events.tsv')
  df = read.table(fName,sep='\t',header=TRUE,row.names=NULL,stringsAsFactors=FALSE)
  df$case_id = rep(x,nrow(df))
  return(df)
},mc.cores=8)
dfl = do.call(rbind,lsv)

#GridSS
gFiles = list.files(paste0(workdir,'Raw Data/Gridss Filtered Events Examples'))
gFiles = as.list(sapply(gFiles,FUN=function(x){strsplit(x,split='_')[[1]][1]},USE.NAMES=FALSE))
gsv = mclapply(gFiles,FUN=function(x){
  fName = paste0(workdir,'Raw Data/Gridss Filtered Events Examples/',x,'_filtered_events.tsv')
  df = read.table(fName,sep='\t',header=TRUE,row.names=NULL,stringsAsFactors=FALSE)
  df$case_id = rep(x,nrow(df))
  return(df)
},mc.cores=8)
dfg = do.call(rbind,gsv)

#Novobreak
nFiles = list.files(paste0(workdir,'Raw Data/Novobreak Filtered Events Examples'))
nFiles = as.list(sapply(nFiles,FUN=function(x){strsplit(x,split='_')[[1]][1]},USE.NAMES=FALSE))
nsv = mclapply(nFiles,FUN=function(x){
  fName = paste0(workdir,'Raw Data/Novobreak Filtered Events Examples/',x,'_filtered_events.tsv')
  df = read.table(fName,sep='\t',header=TRUE,row.names=NULL,stringsAsFactors=FALSE)
  df$case_id = rep(x,nrow(df))
  return(df)
},mc.cores=8)
dfn = do.call(rbind,nsv)

#Manta
mFiles = list.files(paste0(workdir,'Raw Data/Manta Filtered Events Examples'))
mFiles = as.list(sapply(mFiles,FUN=function(x){strsplit(x,split='_')[[1]][1]},USE.NAMES=FALSE))
msv = mclapply(mFiles,FUN=function(x){
  fName = paste0(workdir,'Raw Data/Manta Filtered Events Examples/',x,'_filtered_events.tsv')
  df = read.table(fName,sep='\t',header=TRUE,row.names=NULL,stringsAsFactors=FALSE)
  df$case_id = rep(x,nrow(df))
  return(df)
},mc.cores=8)
dfm = do.call(rbind,msv)

###intersect filtered SV events for each caller with transcripts bed file
mantaSVIntersect = detect_overlaps(df=dfm,bed=bed)
lumpySVIntersect = detect_overlaps(df=dfl,bed=bed)
novobreakSVIntersect = detect_overlaps(df=dfn,bed=bed)
gridssSVIntersect = detect_overlaps(df=dfg,bed=bed)

###save intersected SVs if desired
#write.table(mantaSVIntersect,'',sep='\t',quote=FALSE,col.names=TRUE,row.names=FALSE)
#write.table(lumpySVIntersect,'',sep='\t',quote=FALSE,col.names=TRUE,row.names=FALSE)
#write.table(novobreakSVIntersect,'',sep='\t',quote=FALSE,col.names=TRUE,row.names=FALSE)
#write.table(gridssSVIntersect,'',sep='\t',quote=FALSE,col.names=TRUE,row.names=FALSE)
