##################################################################################################################
### Description: Cluster tumors by RNA singature profile similarity
###              Visualize with heatmaps

library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stats)
workdir = '/rocker-build/gmkf_nbl_somatic/'

### tumor data file
clinData = read.table(paste0(workdir,'Data/tumor_clinical_data.tsv'),sep='\t',stringsAsFactors = FALSE,header=TRUE)

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


### perform heirarchical clustering for gene signatures across all tumors
###cluster signatures together
#remove mes/adrn data for now
pdCluster = pDF %>% filter(!(feature %in% c('Mesenchymal','Slight Mesenchymal','Slight Adrenergic','Adrenergic')))
#get unique gene sets and add to data
rSigs = unique(sapply(as.character(pdCluster$feature),FUN=function(x){strsplit(x,split='[.]')[[1]][1]},USE.NAMES=FALSE))
pdCluster$rsig = sapply(as.character(pdCluster$feature),FUN=function(x){strsplit(x,split='[.]')[[1]][1]},USE.NAMES=FALSE)

#convert bin to score for each sample/gene set
pdClusterdf = data.frame()
for (c in unique(pdCluster$case_id)) {
  fVals = c(c)
  for (r in rSigs) {
    f = pdCluster %>% filter(case_id == c & rsig == r) %>% pull(feature) %>% as.character()
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
  fvdf = t(as.data.frame(fVals,col.names=NULL,row.names=c('case_id',rSigs)))
  pdClusterdf = rbind(pdClusterdf,fvdf)
  rownames(pdClusterdf) = NULL
}

#add sample as row name and transpose
rownames(pdClusterdf) = pdClusterdf$case_id
pdClusterdf = pdClusterdf[,(2:ncol(pdClusterdf))]
pdClusterdf = t(pdClusterdf) %>% as.data.frame()

#hcluster rsig data to get levels order
clustersig = mutate_all(pdClusterdf, function(x) {as.numeric(as.character(x))})
rownames(clustersig) = rownames(pdClusterdf)
clustersig = stats::hclust(dist(clustersig))
clusterDend = plot(clustersig,hang=-0.01,xlab='',main=NULL,sub='',axes=FALSE,ylab=NULL)

#get order of labels
rsigOrder = clustersig$labels[clustersig$order]

#get all gene signatures
gSigs = rsigOrder

#add y-value for plotting
pDF$val = rep(1,nrow(pDF))

#remove samples with no RNA data
pDF = pDF %>% filter(!(case_id %in% swoRNA))

####set factor order
unGroups = unique(pDF$feature)
#create a factor only if feature is present in data
#vector of all possible feature combinations
pref = c(gSigs)
suff = c(apply(expand.grid('NEG',20:1),MARGIN=1,FUN=function(x){paste0(x['Var1'],'.',as.character(as.integer(x['Var2'])))}),'ZERO',
         apply(expand.grid('POS',1:20),MARGIN=1,FUN=function(x){paste0(x['Var1'],'.',as.character(as.integer(x['Var2'])))}),'OUTLIER','NODATA'
)
allFeatureCombo = apply(expand.grid(pref,suff),MARGIN=1,FUN=function(x){paste0(x['Var1'],'.',x['Var2'])})

#check if each feature exists in data, if so, add to vector
#must go in order of pref,suff
gVec = c()
for (g in pref) {
  for (gr in suff) {
    f = paste0(g,'.',gr)
    if (f %in% unGroups) {
      gVec = c(gVec,f)
    }
  }
}

#add in mes/adrn status to feature order
fOrder = c(c('Mesenchymal','Slight Mesenchymal','Slight Adrenergic','Adrenergic'),gVec)

#length(fOrder)
pDF$feature = factor(pDF$feature,levels=fOrder)

#set color schemes into vector
reds <- brewer.pal(9, "Reds") 
reds <- colorRampPalette(reds)(20)
greens <- brewer.pal(9, "Greens")
greens <- colorRampPalette(greens)(20)

pal = c('#E4BF47','#F3DD93','#EEADFF','#D424FF') #colors for mes-adrn status
#in order of feature factor level, assign correct coloring
for (i in 5:length(fOrder)) {
  f = fOrder[i]
  #split feature string
  fSplit = strsplit(f,split='[.]')[[1]]
  #check if outlier or zero
  if (length(fSplit) < 3) {
    #if outlier
    if (fSplit[2] == 'ZERO') {
      pal = c(pal,'white')
    } else if (fSplit[2] == 'OUTLIER') { #if zero
      pal = c(pal,'cyan')
    } else if (fSplit[2] == 'NODATA') { #if no data
      pal = c(pal,'gray40')
    } 
  } else if (length(fSplit) == 3) {
    #get pos or neg direction 
    direction = fSplit[2]
    #get intensity
    level = as.integer(fSplit[3])
    if (direction == 'POS') {
      pal = c(pal,greens[level])
    } else if (direction == 'NEG') {
      pal = c(pal,reds[level])
    }
  }
}

#ensure pallete length matches number of features
length(pal)
length(fOrder)

#set xvalue for each sample
xDict = setNames(clinData$xcord,nm=clinData$case_id)
pDF$xcord = sapply(pDF$case_id,FUN=function(x){xDict[[x]]})

#reset x-coordinate now that no RNA samples are missing
#set each x-coord, going from min value to max value
numSamp = length(unique(pDF$case_id))
minValExclude = c()
for (i in 1:numSamp) {
  #find minimum coordinate and assign i as value
  minVal = min(pDF$xcord[!(pDF$xcord %in% minValExclude)])
  pDF$xcord[pDF$xcord == minVal] = i
  minValExclude = c(minValExclude,i)
}

#get num samples in each mycn/risk grouping
uSample = unique(pDF$case_id)
numGroups = table(tmb[tmb$case_id %in% uSample,]$group2)

#plot all tumors
allSampPlotNoCluster <- ggplot(data=pDF,aes(x=xcord,y=val,fill=feature)) + geom_bar(stat='identity',width=1) + 
  theme_minimal() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank()) + 
  theme(axis.text.y=element_text(size=13),axis.title.y=element_text(size=18,margin=margin(t=0,r=100,b=0,l=0),angle=0)) + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  scale_fill_manual(values=pal) +
  #scale_fill_manual(values=pal,labels = c('Mesenchymal','Slight Mesenchymal','Slight Adrenergic','Adrenergic','No Data',rep('Negative',5),'Zero',rep('Positive',5),'Strong Outlier','No Data')) +
  geom_hline(yintercept=c(0:28)) +
  theme(axis.line=element_blank()) + ylab('RNA Expression Signatures') + 
  theme(legend.margin=margin(t=0,r=0,b=200,l=125)) + 
  geom_vline(xintercept=c(0.5,22.5,53.5,102.5,210.5,201.5)) + 
  coord_cartesian(xlim=c(0.5,201.5),ylim=c(0,28),expand=F,clip="off") + 
  scale_y_continuous(breaks=seq(0.5,27.5,by=1),labels=c(rev(gSigs),'MES/ADRN')) + 
  guides(fill=guide_legend(ncol=1,byrow=TRUE)) +
  theme(legend.position='None')
allSampPlotNoCluster


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
    print(g)
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
    print(ff)
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

# get factors in group and order
unGroupsLR = unique(as.character(pdfLR$feature))
unGroupsIR = unique(as.character(pdfIR$feature))
unGroupsHRN = unique(as.character(pdfHRN$feature))
unGroupsHRM = unique(as.character(pdfHRM$feature))

# create a factor only if feature is present in data
#vector of all possible feature combinations with custom order
customGSigOrder = c('myc','dna_repair',"g2m",'pten_ko','atm_ko','hypoxia','p53','ccnd1','kras','egfr','raf1','ezh2_ko','akt',
                     'erbb2','mek','nfkb','angiogenesis','emt','mtor','pi3k','alk','yap1','wnt','notch','tgfb','hh','gli1')
pref = c(customGSigOrder)
suff = c(apply(expand.grid('NEG',20:1),MARGIN=1,FUN=function(x){paste0(x['Var1'],'.',as.character(as.integer(x['Var2'])))}),'ZERO',
         apply(expand.grid('POS',1:20),MARGIN=1,FUN=function(x){paste0(x['Var1'],'.',as.character(as.integer(x['Var2'])))}),'OUTLIER','NODATA'
)
allFeatureCombo = apply(expand.grid(pref,suff),MARGIN=1,FUN=function(x){paste0(x['Var1'],'.',x['Var2'])})

#check if each feature exists in data, if so, add to vector
#must go in order of pref,suff
gVecAll = lapply(list(unGroupsLR,unGroupsIR,unGroupsHRN,unGroupsHRM),FUN=function(x){
  gVec = c() 
  for (g in pref) {
    for (gr in suff) {
      f = paste0(g,'.',gr)
      if (f %in% x) {
        gVec = c(gVec,f)
      }
    }
  }
  return(gVec)
})
gVecLR = gVecAll[[1]]
gVecIR = gVecAll[[2]]
gVecHRN = gVecAll[[3]]
gVecHRM = gVecAll[[4]]

#check to see which mes/adrn statuses are in each group
c('Mesenchymal','Slight Mesenchymal','Slight Adrenergic','Adrenergic') %in% unGroupsLR
c('Mesenchymal','Slight Mesenchymal','Slight Adrenergic','Adrenergic') %in% unGroupsIR
c('Mesenchymal','Slight Mesenchymal','Slight Adrenergic','Adrenergic') %in% unGroupsHRN
c('Mesenchymal','Slight Mesenchymal','Slight Adrenergic','Adrenergic') %in% unGroupsHRM
#add these features only if they were present
fOrderLR = c(c('Mesenchymal','Slight Mesenchymal','Slight Adrenergic','Adrenergic'),gVecLR)
fOrderIR = c(c('Mesenchymal','Slight Mesenchymal','Slight Adrenergic','Adrenergic'),gVecIR)
fOrderHRN = c(c('Mesenchymal','Slight Mesenchymal','Slight Adrenergic','Adrenergic'),gVecHRN)
fOrderHRM = c(c('Mesenchymal','Slight Adrenergic','Adrenergic'),gVecHRM)

#set features as factor
pdfLR$feature = factor(pdfLR$feature,levels=fOrderLR)
pdfIR$feature = factor(pdfIR$feature,levels=fOrderIR)
pdfHRN$feature = factor(pdfHRN$feature,levels=fOrderHRN)
pdfHRM$feature = factor(pdfHRM$feature,levels=fOrderHRM)

#set color schemes into vector
reds <- brewer.pal(9, "Reds") 
reds <- colorRampPalette(reds)(20)
greens <- brewer.pal(9, "Greens")
greens <- colorRampPalette(greens)(20)

#colors for mes-adrn status
palLR = c('#E4BF47','#F3DD93','#EEADFF','#D424FF') 
palIR = c('#E4BF47','#F3DD93','#EEADFF','#D424FF') 
palHRN = c('#E4BF47','#F3DD93','#EEADFF','#D424FF') 
palHRM = c('#E4BF47','#EEADFF','#D424FF') 

#in order of feature factor level, assign correct coloring
palsAll = lapply(list(fOrderLR,fOrderIR,fOrderHRN),FUN=function(x){
  
  fpal = c()
  
  for (i in 5:length(x)) {
    f = x[i]
    #split feature string
    fSplit = strsplit(f,split='[.]')[[1]]
    #check if outlier or zero
    if (length(fSplit) < 3) {
      #if outlier
      if (fSplit[2] == 'ZERO') {
        fpal = c(fpal,'white')
      } else if (fSplit[2] == 'OUTLIER') { #if zero
        fpal = c(fpal,'cyan')
      } else if (fSplit[2] == 'NODATA') { #if no data
        fpal = c(fpal,'gray40')
      } 
    } else if (length(fSplit) == 3) {
      #get pos or neg direction 
      direction = fSplit[2]
      #get intensity
      level = as.integer(fSplit[3])
      if (direction == 'POS') {
        fpal = c(fpal,greens[level])
      } else if (direction == 'NEG') {
        fpal = c(fpal,reds[level])
      }
    }
  }
  return(fpal)
})
palLR = c(palLR,palsAll[[1]])
palIR = c(palIR,palsAll[[2]])
palHRN = c(palHRN,palsAll[[3]])

#ensure all pallete lengths match number of features
length(palLR) == length(fOrderLR)
length(palIR) == length(fOrderIR)
length(palHRN) == length(fOrderHRN)

#Same process for HRM group done separately
fpal = c()
for (i in 4:length(fOrderHRM)) {
  f = fOrderHRM[i]
  #split feature string
  fSplit = strsplit(f,split='[.]')[[1]]
  #check if outlier or zero
  if (length(fSplit) < 3) {
    #if outlier
    if (fSplit[2] == 'ZERO') {
      fpal = c(fpal,'white')
    } else if (fSplit[2] == 'OUTLIER') { #if zero
      fpal = c(fpal,'cyan')
    } else if (fSplit[2] == 'NODATA') { #if no data
      fpal = c(fpal,'gray40')
    } 
  } else if (length(fSplit) == 3) {
    #get pos or neg direction 
    direction = fSplit[2]
    #get intensity
    level = as.integer(fSplit[3])
    if (direction == 'POS') {
      fpal = c(fpal,greens[level])
    } else if (direction == 'NEG') {
      fpal = c(fpal,reds[level])
    }
  }
}
palHRM = c(palHRM,fpal)
length(palHRM) == length(fOrderHRM)

##### plot
#set xvalue for each sample
xDictLR = setNames(c(1:length(lrOrder)),nm=lrOrder)
xDictIR = setNames(c(1:length(irOrder)),nm=irOrder)
xDictHRN = setNames(c(1:length(hrnOrder)),nm=hrnOrder)
xDictHRM = setNames(c(1:length(hrmOrder)),nm=hrmOrder)
pdfLR$xcord = sapply(pdfLR$case_id,FUN=function(x){xDictLR[[x]]})
pdfIR$xcord = sapply(pdfIR$case_id,FUN=function(x){xDictIR[[x]]})
pdfHRN$xcord = sapply(pdfHRN$case_id,FUN=function(x){xDictHRN[[x]]})
pdfHRM$xcord = sapply(pdfHRM$case_id,FUN=function(x){xDictHRM[[x]]})

lrPlot <- ggplot(data=pdfLR,aes(x=xcord,y=val,fill=feature)) + geom_bar(stat='identity',width=1) + 
  theme_minimal() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank()) + 
  theme(axis.text.y=element_text(size=13),axis.title.y=element_text(size=18,margin=margin(t=0,r=100,b=0,l=0),angle=0)) + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  scale_fill_manual(values=palLR) +
  #scale_fill_manual(values=pal,labels = c('Mesenchymal','Slight Mesenchymal','Slight Adrenergic','Adrenergic','No Data',rep('Negative',5),'Zero',rep('Positive',5),'Strong Outlier','No Data')) +
  geom_hline(yintercept=c(0:28)) +
  theme(axis.line=element_blank()) + ylab('RNA Expression Signatures') + 
  theme(legend.margin=margin(t=0,r=0,b=200,l=125)) + 
  coord_cartesian(xlim=c(0.5,length(lrOrder)+0.5),ylim=c(0,28),expand=F,clip="off") + 
  scale_y_continuous(breaks=seq(0.5,27.5,by=1),labels=c(rev(customGSigOrder),'MES/ADRN')) + 
  guides(fill=guide_legend(ncol=1,byrow=TRUE)) +
  theme(legend.position='None') + geom_vline(xintercept=c(0.5,length(lrOrder)+0.5))
lrPlot

irPlot <- ggplot(data=pdfIR,aes(x=xcord,y=val,fill=feature)) + geom_bar(stat='identity',width=1) + 
  theme_minimal() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank()) + 
  theme(axis.text.y=element_text(size=13),axis.title.y=element_text(size=18,margin=margin(t=0,r=100,b=0,l=0),angle=0)) + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  scale_fill_manual(values=palIR) +
  #scale_fill_manual(values=pal,labels = c('Mesenchymal','Slight Mesenchymal','Slight Adrenergic','Adrenergic','No Data',rep('Negative',5),'Zero',rep('Positive',5),'Strong Outlier','No Data')) +
  geom_hline(yintercept=c(0:28)) +
  theme(axis.line=element_blank()) + ylab('RNA Expression Signatures') + 
  theme(legend.margin=margin(t=0,r=0,b=200,l=125)) + 
  geom_vline(xintercept=c(0.5,length(irOrder)+0.5)) + 
  coord_cartesian(xlim=c(0.5,length(irOrder)+0.5),ylim=c(0,28),expand=F,clip="off") + 
  scale_y_continuous(breaks=seq(0.5,27.5,by=1),labels=c(rev(customGSigOrder),'MES/ADRN')) + 
  guides(fill=guide_legend(ncol=1,byrow=TRUE)) +
  theme(legend.position='None')
irPlot

hrnPlot <- ggplot(data=pdfHRN,aes(x=xcord,y=val,fill=feature)) + geom_bar(stat='identity',width=1) + 
  theme_minimal() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank()) + 
  theme(axis.text.y=element_text(size=13),axis.title.y=element_text(size=18,margin=margin(t=0,r=100,b=0,l=0),angle=0)) + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  scale_fill_manual(values=palHRN) +
  #scale_fill_manual(values=pal,labels = c('Mesenchymal','Slight Mesenchymal','Slight Adrenergic','Adrenergic','No Data',rep('Negative',5),'Zero',rep('Positive',5),'Strong Outlier','No Data')) +
  geom_hline(yintercept=c(0:28)) +
  theme(axis.line=element_blank()) + ylab('RNA Expression Signatures') + 
  theme(legend.margin=margin(t=0,r=0,b=200,l=125)) + 
  geom_vline(xintercept=c(0.5,length(hrnOrder)+0.5)) + 
  coord_cartesian(xlim=c(0.5,length(hrnOrder)+0.5),ylim=c(0,28),expand=F,clip="off") + 
  scale_y_continuous(breaks=seq(0.5,27.5,by=1),labels=c(rev(customGSigOrder),'MES/ADRN')) + 
  guides(fill=guide_legend(ncol=1,byrow=TRUE)) +
  theme(legend.position='None')
hrnPlot

hrmPlot <- ggplot(data=pdfHRM,aes(x=xcord,y=val,fill=feature)) + geom_bar(stat='identity',width=1) + 
  theme_minimal() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank()) + 
  theme(axis.text.y=element_text(size=13),axis.title.y=element_text(size=18,margin=margin(t=0,r=100,b=0,l=0),angle=0)) + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  scale_fill_manual(values=palHRM) +
  #scale_fill_manual(values=pal,labels = c('Mesenchymal','Slight Mesenchymal','Slight Adrenergic','Adrenergic','No Data',rep('Negative',5),'Zero',rep('Positive',5),'Strong Outlier','No Data')) +
  geom_hline(yintercept=c(0:28)) +
  theme(axis.line=element_blank()) + ylab('RNA Expression Signatures') + 
  theme(legend.margin=margin(t=0,r=0,b=200,l=125)) + 
  geom_vline(xintercept=c(0.5,length(hrmOrder)+0.5)) + 
  coord_cartesian(xlim=c(0.5,length(hrmOrder)+0.5),ylim=c(0,28),expand=F,clip="off") + 
  scale_y_continuous(breaks=seq(0.5,27.5,by=1),labels=c(rev(customGSigOrder),'MES/ADRN')) + 
  guides(fill=guide_legend(ncol=1,byrow=TRUE)) +
  theme(legend.position='None')
hrmPlot


