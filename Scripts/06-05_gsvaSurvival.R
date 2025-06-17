##################################################################################################################
### Description: Cluster tumors by RNA singature profile similarity
###              Group by cluster similariy and perform survival analysis

library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stats)
library(survival)
library(survminer)
library(lubridate)
workdir = '/rocker-build/gmkf_nbl_somatic/'

#tumor data file
clinData = read.table(paste0(workdir,'Data/tumor_clinical_data.tsv'),sep='\t',stringsAsFactors = FALSE,header=TRUE)

#import patient clinical data
patCD = read.table(paste0(workdir,'Data/cognbl_gmkf_clindata.csv'),sep=',',stringsAsFactors = FALSE,header=TRUE)
#filter for somatic data
s = clinData$case_id
patCD = patCD %>% filter(usi %in% s)
rownames(patCD) = patCD$usi
patCD = patCD[s,]
#change outcomes to binary values
patCD$vital_status = sapply(patCD$vital_status,FUN=function(x){ifelse(x=='alive',0,1)},USE.NAMES=FALSE)
table(patCD$vital_status)
patCD$vital_status_efs = sapply(1:nrow(patCD),FUN=function(i){ifelse(patCD[i,'event_free_survival_time'] == patCD[i,'overall_survival_time'],
                                                                     0,1)})

#import mes/adrn scores
maDF = read.table(paste0(workdir,'Data/mes_adrn_scores.tsv'),sep='\t',stringsAsFactors=FALSE,header=TRUE,row.names=NULL)
maDF = maDF[,c('SAMPLE','CLASS')]
colnames(maDF) = c('case_id','feature')

#import gsva quantile binning results
#import gsva quantile results
gDF = read.table(paste0(workdir,'Data/gsva_results_quantiles_plot_allgenes.tsv'),sep='\t',
                 stringsAsFactors=FALSE,header=TRUE,row.names=NULL)

#join two tables
pDF = rbind(maDF,gDF[,c('case_id','feature')])

#remove samples not in clinData table
allSamp = unique(pDF$case_id)
toRemove = allSamp[!(allSamp %in% clinData$case_id)]
pDF <- pDF %>% filter(!(case_id %in% toRemove))

#get samples with no RNA data
swRNA = unique(pDF$case_id)
swoRNA = clinData$case_id[!(clinData$case_id %in% swRNA)]

##### perform heirarchical clustering by tumor sig profile using set order of gene signatures for plotting
### plot by risk and MYCN status
#add clinical group to data
pDF$group = sapply(pDF$case_id,FUN=function(x){clinData$group2[clinData$case_id==x]},USE.NAMES=FALSE)

plotGroup = c('Low Risk,Non-Amplified','MYCN-Amplified')
pdfLR = pDF %>% filter(group %in% plotGroup)
plotGroup = c('Intermediate Risk')
pdfIR = pDF %>% filter(group %in% plotGroup)
plotGroup = c('High Risk,Non-Amplified')
pdfHRN = pDF %>% filter(group %in% plotGroup)
plotGroup = c('High Risk,MYCN-amplified')
pdfHRM = pDF %>% filter(group %in% plotGroup)

#change features to numerical values and create new columns for each feature
pdfLRfeatData = setNames(data.frame(matrix(ncol=length(unique(gDF$gene_name))+1,nrow=0)),nm=c('case_id',unique(gDF$gene_name)))
pdfIRfeatData = setNames(data.frame(matrix(ncol=length(unique(gDF$gene_name))+1,nrow=0)),nm=c('case_id',unique(gDF$gene_name)))
pdfHRNfeatData = setNames(data.frame(matrix(ncol=length(unique(gDF$gene_name))+1,nrow=0)),nm=c('case_id',unique(gDF$gene_name)))
pdfHRMfeatData = setNames(data.frame(matrix(ncol=length(unique(gDF$gene_name))+1,nrow=0)),nm=c('case_id',unique(gDF$gene_name)))

#get features
geneFeatures = unique(gDF$gene_name)

#create feature table for all clinical groups
for (c in unique(pdfLR$case_id)) {
  fVals = c(c)
  
  for (g in geneFeatures) {
    f = gDF %>% filter(case_id == c & gene_name == g) %>% pull(feature)
    f1 = strsplit(f,split='[.]')[[1]]
    if (f1[2] == 'ZERO') {
      ff = 0
    } else if (f1[2] == 'OUTLIER') {
      ff = 25
    } else {
      fval = as.numeric(f1[3])
      ff = ifelse(f1[2]=='POS',fval,fval*(-1))
    }
    fVals = c(fVals,ff)
  }
  fvdf = t(as.data.frame(fVals,col.names=NULL,row.names=c('case_id',geneFeatures)))
  pdfLRfeatData = rbind(pdfLRfeatData,fvdf)
  rownames(pdfLRfeatData) = NULL
}
for (c in unique(pdfIR$case_id)) {
  fVals = c(c)
  for (g in geneFeatures) {
    f = gDF %>% filter(case_id == c & gene_name == g) %>% pull(feature)
    f1 = strsplit(f,split='[.]')[[1]]
    if (f1[2] == 'ZERO') {
      ff = 0
    } else if (f1[2] == 'OUTLIER') {
      ff = 25
    } else {
      fval = as.numeric(f1[3])
      ff = ifelse(f1[2]=='POS',fval,fval*(-1))
    }
    fVals = c(fVals,ff)
  }
  fvdf = t(as.data.frame(fVals,col.names=NULL,row.names=c('case_id',geneFeatures)))
  pdfIRfeatData = rbind(pdfIRfeatData,fvdf)
  rownames(pdfIRfeatData) = NULL
}
for (c in unique(pdfHRN$case_id)) {
  fVals = c(c)
  for (g in geneFeatures) {
    f = gDF %>% filter(case_id == c & gene_name == g) %>% pull(feature)
    f1 = strsplit(f,split='[.]')[[1]]
    if (f1[2] == 'ZERO') {
      ff = 0
    } else if (f1[2] == 'OUTLIER') {
      ff = 25
    } else {
      fval = as.numeric(f1[3])
      ff = ifelse(f1[2]=='POS',fval,fval*(-1))
    }
    fVals = c(fVals,ff)
  }
  fvdf = t(as.data.frame(fVals,col.names=NULL,row.names=c('case_id',geneFeatures)))
  pdfHRNfeatData = rbind(pdfHRNfeatData,fvdf)
  rownames(pdfHRNfeatData) = NULL
}
for (c in unique(pdfHRM$case_id)) {
  fVals = c(c)
  for (g in geneFeatures) {
    f = gDF %>% filter(case_id == c & gene_name == g) %>% pull(feature)
    f1 = strsplit(f,split='[.]')[[1]]
    if (f1[2] == 'ZERO') {
      ff = 0
    } else if (f1[2] == 'OUTLIER') {
      ff = 25
    } else {
      fval = as.numeric(f1[3])
      ff = ifelse(f1[2]=='POS',fval,fval*(-1))
    }
    fVals = c(fVals,ff)
  }
  fvdf = t(as.data.frame(fVals,col.names=NULL,row.names=c('case_id',geneFeatures)))
  pdfHRMfeatData = rbind(pdfHRMfeatData,fvdf)
  rownames(pdfHRMfeatData) = NULL
}

## perform heirarchical clustering for each clinical group
#make sure data is numeric
clustLR = mutate_all(pdfLRfeatData[,2:ncol(pdfLRfeatData)], function(x) {as.numeric(as.character(x))})
clustIR = mutate_all(pdfIRfeatData[,2:ncol(pdfIRfeatData)], function(x) {as.numeric(as.character(x))})
clustHRN = mutate_all(pdfHRNfeatData[,2:ncol(pdfHRNfeatData)], function(x) {as.numeric(as.character(x))})
clustHRM = mutate_all(pdfHRMfeatData[,2:ncol(pdfHRMfeatData)], function(x) {as.numeric(as.character(x))})
rownames(clustLR) = pdfLRfeatData$case_id
rownames(clustIR) = pdfIRfeatData$case_id
rownames(clustHRN) = pdfHRNfeatData$case_id
rownames(clustHRM) = pdfHRMfeatData$case_id
# perform clustering
clustersLR = stats::hclust(dist(clustLR))
clustersIR = stats::hclust(dist(clustIR))
clustersHRN = stats::hclust(dist(clustHRN))
clustersHRM = stats::hclust(dist(clustHRM))
lrDend = plot(clustersLR,hang=-0.01,xlab='',main=NULL,sub='',axes=FALSE,ylab=NULL)
irDend = plot(clustersIR,hang=-0.01,xlab='',main=NULL,sub='',axes=FALSE,ylab=NULL)
hrnDend = plot(clustersHRN,hang=-0.01,xlab='',main=NULL,sub='',axes=FALSE,ylab=NULL)
hrmDend = plot(clustersHRM,hang=-0.01,xlab='',main=NULL,sub='',axes=FALSE,ylab=NULL)

### plot RNA heatmaps in order of clusters
#get order of labels
lrOrder = clustersLR$labels[clustersLR$order]
irOrder = clustersIR$labels[clustersIR$order]
hrnOrder = clustersHRN$labels[clustersHRN$order]
hrmOrder = clustersHRM$labels[clustersHRM$order]

####################### Survival by mes-like adrn samples in low risk subgroup
#create groupings
lrg1 = lrOrder[1:32]
lrg2 = lrOrder[33:length(lrOrder)]

#test only samples with clinical data
sData = patCD %>% filter(usi %in% c(lrg1,lrg2))
length(c(lrg1,lrg2))
dim(sData)

#add new column for two groups
lrgDict = c(setNames(rep('MES-like',length(lrg1)),nm=lrg1),setNames(rep('ADRN-like',length(lrg2)),nm=lrg2))
sData$Expression_Group = sapply(sData$usi,FUN=function(x){lrgDict[[x]]},USE.NAMES=FALSE)
sData$vital_status_efs = sapply(1:nrow(sData),FUN=function(i){ifelse(sData[i,'event_free_survival_time'] == sData[i,'overall_survival_time'],
                                                                     0,1)})

#run efs analysis
ep = survfit(Surv(event_free_survival_time,vital_status_efs) ~ Expression_Group, data=sData)
fill = c('Red','Blue')
ggsurvplot(ep,pval=TRUE, risk.table=TRUE, 
           title="Event Free Survival by MES-like status in ADRN low risk patients", palette=fill,
           risk.table.height=.3)
#run os analysis
op = survfit(Surv(overall_survival_time,vital_status) ~ Expression_Group, data=sData)
fill = c('Red','Blue')
ggsurvplot(op,pval=TRUE, risk.table=TRUE, 
           title="Overall Survival Probability by MES-like status in ADRN low risk patients", palette=fill,
           risk.table.height=.3)

####################### Survival by mes-like adrn samples in int risk subgroup
#create groupings
irg1 = irOrder[1:10]
irg2 = irOrder[11:length(irOrder)]

#test only samples with clinical data
sData = patCD %>% filter(usi %in% c(irg1,irg2))
length(c(irg1,irg2))
dim(sData)

#add new column for two groups
irgDict = c(setNames(rep('MES-like ADRN',length(irg1)),nm=irg1),setNames(rep('ADRN',length(irg2)),nm=irg2))
sData$Expression_Group = sapply(sData$usi,FUN=function(x){irgDict[[x]]},USE.NAMES=FALSE)
sData$vital_status_efs = sapply(1:nrow(sData),FUN=function(i){ifelse(sData[i,'event_free_survival_time'] == sData[i,'overall_survival_time'],
                                                                     0,1)})

#run efs analysis
ep = survfit(Surv(event_free_survival_time,vital_status_efs) ~ Expression_Group, data=sData)
fill = c('Red','Blue')
ggsurvplot(ep,pval=TRUE, risk.table=TRUE, 
           title="Event Free Survival by MES-like status in ADRN intermediate risk patients", palette=fill,
           risk.table.height=.3)
#run os analysis
op = survfit(Surv(overall_survival_time,vital_status) ~ Expression_Group, data=sData)
fill = c('Red','Blue')
ggsurvplot(op,pval=TRUE, risk.table=TRUE, 
           title="Overall Survival Probability by MES-like status in ADRN intermediate risk patients", palette=fill,
           risk.table.height=.3)

####################### Survival by mes-like adrn samples in hr-non mycn risk subgroup
#create groupings
hrng1 = hrnOrder[1:6]
hrng2 = hrnOrder[7:19]
hrng3 = hrnOrder[20:25]
hrng4 = hrnOrder[26:length(hrnOrder)]

#test only samples with clinical data
sData = patCD %>% filter(usi %in% c(hrng1,hrng2,hrng3,hrng4))

#add new column for two groups
hrngDict = c(setNames(rep('ADRN-1',length(hrng1)),nm=hrng1),
             setNames(rep('ADRN-2',length(hrng2)),nm=hrng2),
             setNames(rep('MES',length(hrng3)),nm=hrng3),
             setNames(rep('ADRN-3',length(hrng4)),nm=hrng4))
sData$Expression_Group = sapply(sData$usi,FUN=function(x){hrngDict[[x]]},USE.NAMES=FALSE)
sData$vital_status_efs = sapply(1:nrow(sData),FUN=function(i){ifelse(sData[i,'event_free_survival_time'] == sData[i,'overall_survival_time'],
                                                                     0,1)})

#run efs analysis
ep = survfit(Surv(event_free_survival_time,vital_status_efs) ~ Expression_Group, data=sData)
fill =  brewer.pal(4,'Set2')
ggsurvplot(ep,pval=TRUE, risk.table=TRUE, 
           title="Event Free Survival by MES-like status in ADRN high risk non-MYCN amplified patients", palette=fill,
           risk.table.height=.3)
#run os analysis
op = survfit(Surv(overall_survival_time,vital_status) ~ Expression_Group, data=sData)
fill =  brewer.pal(4,'Set2')
ggsurvplot(op,pval=TRUE, risk.table=TRUE, 
           title="Overall Survival Probability by MES-like status in ADRN high risk non-MYCN amplified patients", palette=fill,
           risk.table.height=.3)

####################### Survival by mes-like adrn samples in hr mycn-amp risk subgroup
#create groupings
hrmg1 = hrmOrder[1:3]
hrmg2 = hrmOrder[4:9]
hrmg3 = hrmOrder[10:length(hrmOrder)]

#test only samples with clinical data
sData = patCD %>% filter(usi %in% c(hrmg1,hrmg2,hrmg3))

#add new column for two groups
hrmgDict = c(setNames(rep('ADRN-1',length(hrmg1)),nm=hrmg1),
             setNames(rep('ADRN-2',length(hrmg2)),nm=hrmg2),
             setNames(rep('ADRN-3',length(hrmg3)),nm=hrmg3))
sData$Expression_Group = sapply(sData$usi,FUN=function(x){hrmgDict[[x]]},USE.NAMES=FALSE)
sData$vital_status_efs = sapply(1:nrow(sData),FUN=function(i){ifelse(sData[i,'event_free_survival_time'] == sData[i,'overall_survival_time'],
                                                                     0,1)})

#run efs analysis
ep = survfit(Surv(event_free_survival_time,vital_status_efs) ~ Expression_Group, data=sData)
fill =  brewer.pal(3,'Set2')
ggsurvplot(ep,pval=TRUE, risk.table=TRUE, 
           title="Event Free Survival by MES-like status in ADRN high risk non-MYCN amplified patients", palette=fill,
           risk.table.height=.3)
#run os analysis
op = survfit(Surv(overall_survival_time,vital_status) ~ Expression_Group, data=sData)
fill =  brewer.pal(3,'Set2')
ggsurvplot(op,pval=TRUE, risk.table=TRUE, 
           title="Overall Survival Probability by expression group in high risk non-MYCN amplified patients", palette=fill,
           risk.table.height=.3)


