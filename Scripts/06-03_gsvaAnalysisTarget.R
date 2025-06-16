##################################################################################################################
### Description: Perform gsva analysis using manually curated gene sets related to cancer
###              Gene sets have been acquired from GO, as well as several independent research papers
###              This analysis is for the TARGET neuroblastoma cohort

library(GSVA)
library(dplyr)
library(R.utils)
library(MASS)
library(extremevalues)
library(gplots)
library(RColorBrewer)

workdir = '/rocker-build/gmkf_nbl_somatic/'

#import TARGET clinical data
cdt = read.table(paste0(workdir,'Data/TARGET_harmonized_2018-03-31.csv'),sep=',',
                 header=TRUE,row.names=NULL,stringsAsFactors=FALSE)
colnames(cdt)[1] = 'case_id'

#import gmkf tumor clinical and molecular data
gcd = read.table(paste0(workdir,'Data/tumor_clinical_data.tsv'),sep='\t',stringsAsFactors = FALSE,header=TRUE)

#import TARGET RNA readcount data
expMat = read.table(paste0(workdir,'Data/target_expmat.tsv'),sep='\t',header=TRUE,row.names=1,stringsAsFactors=FALSE)

#keep only samples with RNA data and not in GMKF set
cdt = cdt %>% filter(case_id %in% colnames(expMat))

#keep only HR samples if desired
#hrsamp = cdt$case_id[cdt$COG.Risk.Group == 'High Risk']
#expMat = expMat[,(colnames(expMat) %in% hrsamp)]

#remove duplicates from gabby data
gSamps = gcd$case_id
cdt = cdt %>% filter(!(case_id %in% gSamps))
expMat = expMat[,(!(colnames(expMat) %in% gSamps))]

#list of coding genes
codGenes = readRDS(paste0(workdir,'References/hg38_coding_genes.RDS'))

#filter expression matrix by coding genes only
expMatCod = expMat[rownames(expMat) %in% codGenes,]

#import gene sets
geneSets = setNames(list(scan(paste0(workdir,'Data/gene_signatures/hallmark_kras.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/hallmark_myc_v2.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/hallmark_dna_repair.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/hallmark_angiogenesis.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/hallmark_emt.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/hallmark_g2m.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/hallmark_hh.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/hallmark_hypoxia.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/hallmark_mtor.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/hallmark_nfkb.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/hallmark_notch.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/hallmark_p53.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/hallmark_pi3k.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/hallmark_tgfb.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/hallmark_wnt.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/raf1.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/alk.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/akt.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/erbb2.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/ccnd1.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/egfr.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/gli1.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/mek.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/ezh2_ko.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/pten_ko.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/yap1.txt'),character()),
                         scan(paste0(workdir,'Data/gene_signatures/atm_ko.txt'),character())
),
c('kras','myc','dna_repair','angiogenesis','emt','g2m','hh','hypoxia','mtor','nfkb','notch','p53',
  'pi3k','tgfb','wnt','raf1','alk','akt','erbb2','ccnd1','egfr','gli1','mek','ezh2_ko','pten_ko',
  'yap1','atm_ko'))

## run GSVA
res = gsva(expMatCod %>% as.matrix(), geneSets)

## determine quantiles for plotting, identify outliers in each gene signature
outliersByGene = list()
outliersModelRho = list()
outliersModelFlim = list()

#set parameters for calculating outliers
i = 2
r = c(13,13)
f = c(0.3,0.8)
#calculate outliers
outliers = getOutliers(res[i,],method='I',rho=r,FLim=f,distribution='normal')
print(names(geneSets)[i])
print(outliers$nOut)
#plot outliers
outlierPlot(res[i,], outliers, mode="qq")

#calculate outliers for each gene set
outliersByGene[[names(geneSets)[i]]] = names(outliers$iRight)
outliersModelRho[[names(geneSets)[i]]] = r
outliersModelFlim[[names(geneSets)[i]]] = f

outliersData = setNames(list(outliersByGene,outliersModelRho,outliersModelFlim),
                        nm=c('samples','rho_params','Flim_params'))

## look at outliers with all data at once
resV <- as.vector(res)
r = c(10,10)
f = c(0.1,0.99)
outliers = getOutliers(resV,method='I',rho=r,FLim=f,distribution='normal')
outlierPlot(resV, outliers, mode="qq")


### create a data frame of all samples for plotting
#split all data >0 and <0 into 20 bins
#all genes used at once
#outliers will get own label
#set limits for quantiles

#create a function that takes a vector, sequence of pos/neg bins, and returns output as vector of bins
placeBins <- function(v,pb,nb) {
  #pb: vector of bin limits for positive values (increasing)
  #nb: vector of bin limits for negative values (decreasing)
  
  #v = res_copy_list[[2]]
  v = round(v,2)
  xc = as.character(v)
  
  #set zero bin
  xc[v == 0] = 'ZERO'
  
  #assign negative bins
  for (i in 21:2) {
    xc[v >= nb[i] & v < nb[i-1]] = paste0('NEG','.',as.character(i-1))
  }
  
  #assign positive bins
  for (i in 1:20) {
    xc[v > pb[i] & v <= pb[i+1]] = paste0('POS','.',as.character(i))
  }
  
  return(xc)
}


pos = resV[resV > 0]
neg = resV[resV < 0]
posBin = round(seq(0,max(pos),length.out=21),3)
negBin = rev(round(seq(min(neg),0,length.out=21),3))

#assign each sample/gene a bin in a long format df
res_copy = res
res_copy_list = split(res_copy,1:nrow(res_copy)) #split df to lists of rows
binList = lapply(res_copy_list,FUN=placeBins,pb=posBin,nb=negBin)
#create long format df for plotting
s = colnames(res)
df = data.frame()
for (i in 1:27) {
  gene = rownames(res)[i]
  v = binList[[i]]
  v = sapply(v,FUN=function(x){paste0(gene,'.',x)},USE.NAMES=FALSE)
  for (k in 1:length(s)) {
    df = rbind(df,data.frame(case_id=s[k],feature=v[k],gene_name = gene))
  }
}
df$feature = as.character(df$feature)
df$gene_name = as.character(df$gene_name)
df$case_id = as.character(df$case_id)

#add outliers
for (g in rownames(res)) {
  oSamp = outliersData$samples[[g]]
  #find sample in df and change to format {[[g]]OUTLIER}
  for (s in oSamp)
    df[df$case_id == s & df$gene_name ==g,'feature'] = paste0(g,'.','OUTLIER')
}

### create a data frame of all samples for plotting
#split all data >0 and <0 into 5 quantiles
#outliers will get own label
dfRNA = data.frame()

for (i in c(1:nrow(res))) {
  
  row = round(res[i,],2)
  gene = rownames(res)[i]
  outlierSamp = outliersData[['samples']][[gene]]
  
  pos = row[row > 0]
  pq = quantile(pos, probs=seq(0,1,by=0.2))
  s1qp = names(pos[pos >= pq[1] & pos < pq[2]])
  s2qp = names(pos[pos >= pq[2] & pos < pq[3]])
  s3qp = names(pos[pos >= pq[3] & pos < pq[4]])
  s4qp = names(pos[pos >= pq[4] & pos < pq[5]])
  s5qp = names(pos[pos >= pq[5] & pos <= pq[6]])
  
  neg = row[row < 0]
  nq = quantile(neg, probs=seq(0,1,by=0.2))
  s1qn = names(neg[neg <= nq[6] & neg > nq[5]])
  s2qn = names(neg[neg <= nq[5] & neg > nq[4]])
  s3qn = names(neg[neg <= nq[4] & neg > nq[3]])
  s4qn = names(neg[neg <= nq[3] & neg > nq[2]])
  s5qn = names(neg[neg <= nq[2] & neg >= nq[1]])
  
  zero = names(row[row==0])
  
  df = data.frame()
  
  for (s in s1qp) {
    df = rbind(df,data.frame(case_id=s,quantile='pos1',gene=gene))
  }
  for (s in s2qp) {
    df = rbind(df,data.frame(case_id=s,quantile='pos2',gene=gene))
  }
  for (s in s3qp) {
    df = rbind(df,data.frame(case_id=s,quantile='pos3',gene=gene))
  }
  for (s in s4qp) {
    df = rbind(df,data.frame(case_id=s,quantile='pos4',gene=gene))
  }
  for (s in s5qp) {
    df = rbind(df,data.frame(case_id=s,quantile='pos5',gene=gene))
  }
  for (s in zero) {
    df = rbind(df,data.frame(case_id=s,quantile='zero',gene=gene))
  }
  for (s in s1qn) {
    df = rbind(df,data.frame(case_id=s,quantile='neg1',gene=gene))
  }
  for (s in s2qn) {
    df = rbind(df,data.frame(case_id=s,quantile='neg2',gene=gene))
  }
  for (s in s3qn) {
    df = rbind(df,data.frame(case_id=s,quantile='neg3',gene=gene))
  }
  for (s in s4qn) {
    df = rbind(df,data.frame(case_id=s,quantile='neg4',gene=gene))
  }
  for (s in s5qn) {
    df = rbind(df,data.frame(case_id=s,quantile='neg5',gene=gene))
  }
  
  #change from factor to character
  df$quantile = as.character(df$quantile)
  
  #label outliers
  df[df$case_id %in% outlierSamp,'quantile'] = 'outlier'
  
  #join to full df
  dfRNA = rbind(dfRNA,df)
}

#confirm dimensions of plotting data frame and gsva scores
dfDim = dim(dfRNA)
geneSetSums = apply(res,MARGIN=1,sum)
