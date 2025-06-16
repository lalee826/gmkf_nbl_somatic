##################################################################################################################
### Description: Survival analysis based on Mutation Signature clustering results
###              

library(survival)
library(survminer)
library(lubridate)
library(dplyr)
workdir = '/rocker-build/gmkf_nbl_somatic/'

## import coding mutations
# The following samples were found to have likely tumor/normal mismatches or low NGS checkmate correlation in a separate analysis
excludeSamp = readRDS(c(paste0(workdir,'Data/excluded_samples.RDS')))

### load all coding mut data ###
codMuts = read.table(paste0(workdir,'Data/coding_mutations_annotated.tsv'),sep='\t',header=TRUE,row.names=NULL,stringsAsFactors=FALSE)
#remove excluded samples
codMuts <- codMuts %>% filter(!(Tumor_Sample_Barcode %in% excludeSamp))

## import tumor clinical data table
clinData = read.csv(paste0(workdir,'Data/tumor_clinical_data.tsv'),sep='\t',header=TRUE,row.names=NULL,stringsAsFactors=FALSE)

###### import patient clinical data
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

##################### Survival by mutation signature cluster
#get sig profile and clustering results
csigs = read.table(paste0(workdir,'Data/cosmicv3_sig_contribution.tsv'),sep='\t',stringsAsFactors=FALSE,header=TRUE,row.names=1)
msCluster = read.table(paste0(workdir,'Data/mutsig_cluster_results.tsv'),sep='\t',stringsAsFactors=FALSE,header=TRUE,row.names=NULL)
#make data frame with tumor USI and cluster result
sigClust = sapply(rownames(csigs),FUN=function(x){msCluster %>% filter(sample_id == x) %>% pull(cluster)},USE.NAMES=FALSE)
csigs$cluster = sigClust+1
rownames(csigs) = sapply(rownames(csigs),FUN=function(x){clinData %>% filter(Tumor_Sample_Barcode == x) %>% pull(case_id)},USE.NAMES=FALSE)
csigs = csigs[patCD$usi,]

#make data
sData = cbind(csigs$cluster,patCD) %>% dplyr::rename('MutSig_cluster' = 'csigs$cluster')
sDataMYCN = sData %>% filter(mycn_status == 'amplified')
sDataNonMyc = sData %>% filter(mycn_status != 'amplified')
sDataHR = sData %>% filter(risk == 'high' & mycn_status == 'not amplified')
sDataIRLR = sData %>% filter(risk != 'high')
sDataIR = sData %>% filter(risk == 'intermediate')
sDataLR = sData %>% filter(risk == 'low')

### all data
#run efs analysis
ep = survfit(Surv(event_free_survival_time,vital_status_efs) ~ MutSig_cluster, data=sData)
fill = brewer.pal(length(unique(csigs$cluster)),'Set2')
ggsurvplot(ep,pval=TRUE, risk.table=TRUE, 
           title="Event Free Survival Probability by Mutation Signature Cluster Across all Samples", palette=fill,
           risk.table.height=.3)
#run os analysis
op = survfit(Surv(overall_survival_time,vital_status) ~ MutSig_cluster, data=sData)
fill = brewer.pal(length(unique(csigs$cluster)),'Set2')
ggsurvplot(op,pval=TRUE, risk.table=TRUE, 
           title="Overall Survival Probability by Mutation Signature Cluster Across all Samples", palette=fill,
           risk.table.height=.3)

### MYCN patients
#run efs analysis
ep = survfit(Surv(event_free_survival_time,vital_status_efs) ~ MutSig_cluster, data=sDataMYCN)
fill = brewer.pal(length(unique(csigs$cluster)),'Set2')
ggsurvplot(ep,pval=TRUE, risk.table=TRUE, 
           title="Event Free Survival Probability by Mutation Signature Cluster in MYCN-Positive Patients", palette=fill,
           risk.table.height=.3)
#run os analysis
op = survfit(Surv(overall_survival_time,vital_status) ~ MutSig_cluster, data=sDataMYCN)
fill = brewer.pal(length(unique(csigs$cluster)),'Set2')
ggsurvplot(op,pval=TRUE, risk.table=TRUE, 
           title="Overall Survival Probability by Mutation Signature Cluster in MYCN-Positive Patients", palette=fill,
           risk.table.height=.3)

### NON-MYCN amplified
#run efs analysis
ep = survfit(Surv(event_free_survival_time,vital_status_efs) ~ MutSig_cluster, data=sDataNonMyc)
fill = brewer.pal(length(unique(csigs$cluster)),'Set2')
ggsurvplot(ep,pval=TRUE, risk.table=TRUE, 
           title="Event Free Survival Probability by Mutation Signature Cluster in MYCN-Negative Patients", palette=fill,
           risk.table.height=.3)
#run os analysis
op = survfit(Surv(overall_survival_time,vital_status) ~ MutSig_cluster, data=sDataNonMyc)
fill = brewer.pal(length(unique(csigs$cluster)),'Set2')
ggsurvplot(op,pval=TRUE, risk.table=TRUE, 
           title="Overall Survival Probability by Mutation Signature Cluster in MYCN-Negative Patients", palette=fill,
           risk.table.height=.3)