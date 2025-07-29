##################################################################################################################
### Description: Generate custom oncoplot with all pathogenic events (SNVs, SVs, CNVs)
###      

library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(parallel)
library(data.table)
library(stats)
workdir = '/rocker-build/gmkf_nbl_somatic/'

### import pathogenic events file
allPath = read.table(paste0(workdir,'Data/gmkf_tumor_allPathEvents.tsv'),sep='\t',
                     stringsAsFactors=FALSE,header=TRUE)

#tumor data file
clinData = read.table(paste0(workdir,'Data/tumor_clinical_data.tsv'),sep='\t',stringsAsFactors = FALSE,header=TRUE)



### combine genes into gene families

#combine all DNMT genes
allPath$Gene = sapply(allPath$Gene,FUN=function(x){ifelse(x %in% c('DNMT3A','DNMT3B'),'DNMT Family',x)},USE.NAMES=FALSE)
#combine all KMT2 genes
allPath$Gene = sapply(allPath$Gene,FUN=function(x){ifelse(x %in% c('KMT2A','KMT2A','KMT2D'),'KMT2 Family',x)},USE.NAMES=FALSE)
#combine all FANC genes
allPath$Gene = sapply(allPath$Gene,FUN=function(x){ifelse(x %in% c('FANCA','FANCB','FANCD2','FANCE','FANCG','FANCI'),'FANC Family',x)},USE.NAMES=FALSE)
#combine all RAS genes
allPath$Gene = sapply(allPath$Gene,FUN=function(x){ifelse(x %in% c('KRAS','NRAS','HRAS'),'RAS Family',x)},USE.NAMES=FALSE)
#combine all PIK32 genes
allPath$Gene = sapply(allPath$Gene,FUN=function(x){ifelse(x %in% c('PIK3C2A','PIK3C2B'),'PIK3C2 Family',x)},USE.NAMES=FALSE)
### combine all SMAD genes
allPath$Gene = sapply(allPath$Gene,FUN=function(x){ifelse(x %in% c('SMAD2','SMAD4'),'SMAD Family',x)},USE.NAMES=FALSE)
#combine all SMARC genes
allPath$Gene = sapply(allPath$Gene,FUN=function(x){ifelse(x %in% c('SMARCB1','SMARCC1'),'SMARC Family',x)},USE.NAMES=FALSE)
#combine all FGFR genes
allPath$Gene = sapply(allPath$Gene,FUN=function(x){ifelse(x %in% c("FGFR1","FGFR2","FGFR4"),'FGFR Family',x)},USE.NAMES=FALSE)
#combine all ATM gene names
allPath$Gene = sapply(allPath$Gene,FUN=function(x){ifelse(x %in% c("ATM","ATM;C11orf65"),'ATM',x)},USE.NAMES=FALSE)
#combine all AXIN genes
allPath$Gene = sapply(allPath$Gene,FUN=function(x){ifelse(x %in% c("AXIN1","AXIN2"),'AXIN Family',x)},USE.NAMES=FALSE)
#combine all BRCA genes
allPath$Gene = sapply(allPath$Gene,FUN=function(x){ifelse(x %in% c("BRCA1","BRCA2"),'BRCA Family',x)},USE.NAMES=FALSE)
#combine all ERCC genes
allPath$Gene = sapply(allPath$Gene,FUN=function(x){ifelse(x %in% c("ERCC2","ERCC3",'ERCC4','ERCC5'),'ERCC Family',x)},USE.NAMES=FALSE)
#combine all FGF genes
allPath$Gene = sapply(allPath$Gene,FUN=function(x){ifelse(x %in% c("FGF14","FGF19",'FGF23'),'FGF Family',x)},USE.NAMES=FALSE)
#combine all GATA genes
allPath$Gene = sapply(allPath$Gene,FUN=function(x){ifelse(x %in% c("GATA2","GATA3"),'GATA Family',x)},USE.NAMES=FALSE)
#combine all POL genes
allPath$Gene = sapply(allPath$Gene,FUN=function(x){ifelse(x %in% c("POLD1","POLE",'POLQ'),'POL Family',x)},USE.NAMES=FALSE)
#combine all PTPN genes
allPath$Gene = sapply(allPath$Gene,FUN=function(x){ifelse(x %in% c("PTPN11","PTPN13",'PTPN14'),'PTPN Family',x)},USE.NAMES=FALSE)
#combine all RAD genes
allPath$Gene = sapply(allPath$Gene,FUN=function(x){ifelse(x %in% c("RAD50","RAD51",'RAD51C','RAD51D','RAD54B','TH2LCRR;RAD50'),'RAD Family',x)},USE.NAMES=FALSE)


### Data wrangling for plotting
#split fusions into multiple rows
allPathFusion = allPath[sapply(allPath$Gene,FUN=function(x){length(strsplit(x,'-')[[1]])==2},USE.NAMES=FALSE),]
allPath = allPath[sapply(allPath$Gene,FUN=function(x){length(strsplit(x,'-')[[1]])==1},USE.NAMES=FALSE),]

for (i in 1:nrow(allPathFusion)) {
  row = allPathFusion[i,]
  gSplit = strsplit(row %>% pull('Gene'),'-')[[1]]
  for (g in gSplit) {
    newRow = row
    newRow['Gene'] = g
    allPath = rbind(allPath,newRow)
  }
}

#get all samples
allSamp = clinData$case_id
#get all genes we will use
allGene = unique(allPath$Gene)

#create new sample/gene combo with new gene groups
allPath$case_gene = sapply(1:nrow(allPath),FUN=function(i){
  return(paste(c(allPath[i,'case_id'],allPath[i,'Gene']),collapse='-'))
},USE.NAMES=FALSE)

#remove duplicated rows
allPath = allPath[!(duplicated(allPath$case_gene)),]


### plot genes with count greater than given value
n = 2
geneCount = sort(table(allPath$Gene),decreasing=TRUE)
sum(geneCount > n)
## filter data for genes w count > n
allGene = names(geneCount)[geneCount > n]
allPathG = allPath %>% filter(Gene %in% allGene)

#add in blank data for all samples/genes if combo does not exist in data
newDF = data.frame()
for (s in allSamp) {
  for (g in allGene) {
    sg = paste0(s,'-',g)
    if (sg %in% allPathG$case_gene) {
      next
    } else {
      newRow = data.frame(case_id=s,Gene=g,Type='ND',case_gene=sg)
      newDF = rbind(newDF,newRow)
    }
  }
}

#add in risk group data
riskDict = setNames(c('Low Risk','Int Risk','High Risk'),nm=c('low risk','intermediate risk','high risk'))
for (s in allSamp) {
  risk = riskDict[[clinData$Risk[clinData$case_id == s]]]
  newRow = data.frame(case_id=s,Gene='Risk Group',Type=risk,case_gene=paste0(s,'-Risk'))
  newDF = rbind(newDF,newRow)
}

#combine all data for plotting
allPathG = rbind(allPathG,newDF)

#change all protein coding ablations to same type 
allPathG$Type = sapply(allPathG$Type,FUN=function(x){
                      ifelse(x %in% c("frameshift_variant","stop_gained","start_lost"),'Frameshift/Protein Coding Loss',x)},
                      USE.NAMES=FALSE)
#change missense_variant to Missense
allPathG$Type = sapply(allPathG$Type,FUN=function(x){
                      ifelse(x == 'missense_variant','Missense',x)},
                      USE.NAMES=FALSE)
#change splice variants to Splice Variant
allPathG$Type = sapply(allPathG$Type,FUN=function(x){
  ifelse(x %in% c("missense_variant,splice_region_variant","splice_acceptor_variant","splice_donor_variant"),
         'Splice Variant',x)},
  USE.NAMES=FALSE)


### set x order and add coordinates to table
#* xorder will be set within each risk group by iteratively going through each top gene for each sample
#* split into 4 groups: HR MYCN+, HR MYCN-, INT RISK, LOW RISK(1 MYCN+ sample first)

#order genes with count >= 2
gOrder = names(geneCount[geneCount > n])
gLen = length(gOrder)

## HR MYCN+
#build binary matrix of gene pos for each sample
hrms = clinData$case_id[clinData$Risk == 'high risk' & clinData$Mycn_Status == 'Amplified']
hrmMat = matrix(nrow=0,ncol=gLen)
for (s in hrms) {
  sDat = allPathG %>% filter(case_id == s & Gene != 'Risk Group')
  hrmMat = rbind(hrmMat,sapply(gOrder,FUN=function(x){
    gDat = sDat %>% filter(Gene == x) %>% pull(Type)
    if (length(gDat) > 1) {
      if (sum(gDat != 'ND') > 0) {
        return(1)
      } else {
        return(0)
      }
    } else {
      return(ifelse(gDat == 'ND',0,1))
    }},USE.NAMES=TRUE))
}
hrmDF = as.data.frame(hrmMat)
rownames(hrmDF) = hrms

#print this string to copy/paste into order function
cat(paste0('-hrmDF$"',gOrder, '"'), sep = ",")
#arrange order by path gene count
hrmDF = hrmDF[order(-hrmDF$"MYCN",-hrmDF$"ALK",-hrmDF$"CIC",-hrmDF$"ATRX",-hrmDF$"SHANK2",-hrmDF$"TERT",
            -hrmDF$"CDKN2A",-hrmDF$"FGF Family",-hrmDF$"CREBBP",-hrmDF$"PTPRD",-hrmDF$"ZFHX3",
            -hrmDF$"CCND1",-hrmDF$"FGFR Family",-hrmDF$"PTPN Family",-hrmDF$"CNTNAP2",-hrmDF$"DDR2",
            -hrmDF$"CHD2",-hrmDF$"DDX10",-hrmDF$"GATA Family",-hrmDF$"PTEN",-hrmDF$"RAS Family",
            -hrmDF$"SKA3",-hrmDF$"CDK12",-hrmDF$"DNMT Family",-hrmDF$"TP53",-hrmDF$"ZNF429"), ]
hrmOrder = rownames(hrmDF)


## HR
#build binary matrix of gene pos for each sample
hrs = clinData$case_id[clinData$Risk == 'high risk' & clinData$Mycn_Status == 'Non-Amplified']
hrMat = matrix(nrow=0,ncol=gLen)
for (s in hrs) {
  sDat = allPathG %>% filter(case_id == s & Gene != 'Risk Group')
  hrMat = rbind(hrMat,sapply(gOrder,FUN=function(x){
    gDat = sDat %>% filter(Gene == x) %>% pull(Type)
    if (length(gDat) > 1) {
      if (sum(gDat != 'ND') > 0) {
        return(1)
      } else {
        return(0)
      }
    } else {
      return(ifelse(gDat == 'ND',0,1))
    }},USE.NAMES=TRUE))
}
hrDF = as.data.frame(hrMat)
rownames(hrDF) = hrs

#print this string to copy/paste into order function
cat(paste0('-hrDF$"',gOrder, '"'), sep = ",")
#arrange order by path gene count
hrDF = hrDF[order(-hrDF$"MYCN",-hrDF$"ALK",-hrDF$"CIC",-hrDF$"ATRX",-hrDF$"SHANK2",-hrDF$"TERT",-hrDF$"CDKN2A",
                  -hrDF$"FGF Family",-hrDF$"CREBBP",-hrDF$"PTPRD",-hrDF$"ZFHX3",-hrDF$"CCND1",-hrDF$"FGFR Family",
                  -hrDF$"PTPN Family",-hrDF$"CNTNAP2",-hrDF$"DDR2",-hrDF$"CHD2",-hrDF$"DDX10",-hrDF$"GATA Family",
                  -hrDF$"PTEN",-hrDF$"RAS Family",-hrDF$"SKA3",-hrDF$"CDK12",-hrDF$"DNMT Family",
                  -hrDF$"TP53",-hrDF$"ZNF429"), ]
hrOrder = rownames(hrDF)


## IR
#build binary matrix of gene pos for each sample
irs = clinData$case_id[clinData$Risk == 'intermediate risk']
irMat = matrix(nrow=0,ncol=gLen)
for (s in irs) {
  sDat = allPathG %>% filter(case_id == s & Gene != 'Risk Group')
  irMat = rbind(irMat,sapply(gOrder,FUN=function(x){
    gDat = sDat %>% filter(Gene == x) %>% pull(Type)
    if (length(gDat) > 1) {
      if (sum(gDat != 'ND') > 0) {
        return(1)
      } else {
        return(0)
      }
    } else {
      return(ifelse(gDat == 'ND',0,1))
    }},USE.NAMES=TRUE))
}
irDF = as.data.frame(irMat)
rownames(irDF) = irs

#print this string to copy/paste into order function
cat(paste0('-irDF$"',gOrder, '"'), sep = ",")
#arrange order by path gene count
irDF = irDF[order(-irDF$"MYCN",-irDF$"ALK",-irDF$"CIC",-irDF$"ATRX",-irDF$"SHANK2",-irDF$"TERT",
                  -irDF$"CDKN2A",-irDF$"FGF Family",-irDF$"CREBBP",-irDF$"PTPRD",-irDF$"ZFHX3",
                  -irDF$"CCND1",-irDF$"FGFR Family",-irDF$"PTPN Family",-irDF$"CNTNAP2",-irDF$"DDR2",
                  -irDF$"CHD2",-irDF$"DDX10",-irDF$"GATA Family",-irDF$"PTEN",-irDF$"RAS Family",
                  -irDF$"SKA3",-irDF$"CDK12",-irDF$"DNMT Family",-irDF$"TP53",-irDF$"ZNF429"), ]
irOrder = rownames(irDF)


## LR
#build binary matrix of gene pos for each sample
lrs = clinData$case_id[clinData$Risk == 'low risk']
lrMat = matrix(nrow=0,ncol=gLen)
for (s in lrs) {
  sDat = allPathG %>% filter(case_id == s & Gene != 'Risk Group')
  lrMat = rbind(lrMat,sapply(gOrder,FUN=function(x){
    gDat = sDat %>% filter(Gene == x) %>% pull(Type)
    if (length(gDat) > 1) {
      if (sum(gDat != 'ND') > 0) {
        return(1)
      } else {
        return(0)
      }
    } else {
      return(ifelse(gDat == 'ND',0,1))
    }},USE.NAMES=TRUE))
}
lrDF = as.data.frame(lrMat)
rownames(lrDF) = lrs

#print this string to copy/paste into order function
cat(paste0('-lrDF$"',gOrder, '"'), sep = ",")
#arrange order by path gene count
lrDF = lrDF[order(-lrDF$"MYCN",-lrDF$"ALK",-lrDF$"CIC",-lrDF$"ATRX",-lrDF$"SHANK2",-lrDF$"TERT",
                  -lrDF$"CDKN2A",-lrDF$"FGF Family",-lrDF$"CREBBP",-lrDF$"PTPRD",-lrDF$"ZFHX3",
                  -lrDF$"CCND1",-lrDF$"FGFR Family",-lrDF$"PTPN Family",-lrDF$"CNTNAP2",-lrDF$"DDR2",
                  -lrDF$"CHD2",-lrDF$"DDX10",-lrDF$"GATA Family",-lrDF$"PTEN",-lrDF$"RAS Family",
                  -lrDF$"SKA3",-lrDF$"CDK12",-lrDF$"DNMT Family",-lrDF$"TP53",-lrDF$"ZNF429"), ]
lrOrder = rownames(lrDF)


## get final sample order and check length
sampleOrder = c(hrmOrder,hrOrder,irOrder,lrOrder)
length(sampleOrder)==348

## add x coordinates and fill value to each row in graph data
xCordDict = setNames(c(1:348),nm=sampleOrder)
allPathG$xcord = sapply(allPathG$case_id,FUN=function(x){xCordDict[[x]]},USE.NAMES=FALSE)
allPathG$val = rep(1,nrow(allPathG))

### setting levels and fill colors
#create fill dictionary by mut type
specCols = brewer.pal(11,'Spectral')
darkCols = brewer.pal(8,'Dark2')
rybCols = brewer.pal(8,'RdYlBu')
redCols = brewer.pal(9,'Reds')
purdCols = brewer.pal(9,'PuRd')
grCols = brewer.pal(9,'Greens')
greyCols = brewer.pal(9,'Greys')
yobCols = brewer.pal(9,'YlOrBr')
blCols = brewer.pal(9,'Blues')
purpCols = brewer.pal(9,'Purples')
orCols = brewer.pal(9,'OrRd')
yorCols = brewer.pal(9,'YlOrRd')
s2Cols = brewer.pal(8,'Set2')
pbgCols = brewer.pal(9,'PuBuGn')

fillDict = setNames(c('hotpink1',specCols[3],'gray50','firebrick1','indianred',
                      blCols[8],blCols[6],purpCols[9],purpCols[8],purpCols[7],
                      yorCols[4],pbgCols[7],grCols[6],redCols[6],
                      'gray95',specCols[1],yorCols[5],yorCols[3]),
                    nm = c("Missense","Frameshift/Protein Coding Loss",
                           "Splice Variant","Duplication",
                           "Partial Duplication (Lof)","Deletion","Deletion/Translocation",
                           "Inversion/Translocation","Inversion","Inversion/Deletion",
                           "Insertion","Translocation","Fusion","Amplification",
                           "ND","High Risk","Int Risk","Low Risk"))

## create gene/type combinations and set fill colors
fill = c()
levelsFeat = c()
for (g in c('Risk Group',allGene)) {
  gd = allPathG %>% filter(Gene == g)
  types = unique(gd$Type)
  for (t in types) {
    gt = paste0(g,'-',t)
    fCol = fillDict[[t]]
    levelsFeat = c(levelsFeat,gt)
    fill = c(fill,fCol)
  }
}

## set levels for gene/features
allPathG$Feature = sapply(1:nrow(allPathG),FUN=function(i){paste0(allPathG[i,'Gene'],'-',allPathG[i,'Type'])})
allPathG$Feature = factor(allPathG$Feature,levels=levelsFeat)


#------------------------------------------------------------------------------
### create oncoplot
p <- ggplot(data=allPathG,aes(x=xcord,y=val,fill=Feature)) +
  geom_bar(stat='identity',width=0.90,color='gray90',size=0.4) +
  theme_minimal() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank()) +
  theme(axis.ticks.y=element_blank(),axis.title.y=element_blank()) + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  coord_cartesian(xlim=c(0.5,348.5),ylim=c(0,length(unique(allPathG$Gene))),expand=F,clip="off") + 
  scale_y_continuous(breaks=seq(0.5,length(unique(allPathG$Gene))-0.5,by=1),labels=c(rev(allGene),'Risk Group')) +
  scale_fill_manual(values=fill) + 
  theme(legend.position='None')

grid.newpage()
gt <- ggplotGrob(p)
fig_size <- c(19.5,6) # inches
margin <- unit(2, "line")
gt$vp <- viewport(width = unit(fig_size[1], "in") - margin, height=unit(fig_size[2],"in")- margin)
grid.draw(gt)


#-----------------------------------------------------------------------
#create plot for legend
'''
***********************************************************************
NOTE: plot is incorrect, just using the legend from this figure
***********************************************************************
'''
allDatGLeg = allPathG
allDatGLeg$Type = factor(allDatGLeg$Type,
                         levels = c("Missense","Frameshift/Protein Coding Loss",
                                    "Splice Variant","Duplication",
                                    "Partial Duplication (Lof)","Deletion","Deletion/Translocation",
                                    "Inversion/Translocation","Inversion","Inversion/Deletion",
                                    "Insertion","Translocation","Fusion","Amplification",
                                    "ND","High Risk","Int Risk","Low Risk"))
fillLeg = setNames(c('hotpink1',specCols[3],'gray50','firebrick1','indianred',
            blCols[8],blCols[6],purpCols[9],purpCols[8],purpCols[7],
            yorCols[4],pbgCols[7],grCols[6],redCols[6],
            'gray95',specCols[1],yorCols[5],yorCols[3]),
            nm= c("Missense","Frameshift/Protein Coding Loss",
                  "Splice Variant","Duplication",
                  "Partial Duplication (Lof)","Deletion","Deletion/Translocation",
                  "Inversion/Translocation","Inversion","Inversion/Deletion",
                  "Insertion","Translocation","Fusion","Amplification",
                  "ND","High Risk","Int Risk","Low Risk"))
p <- ggplot(data=allDatGLeg,aes(x=xcord,y=val,fill=Type)) +
  geom_bar(stat='identity',width=1) +
  theme_minimal() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank()) +
  theme(axis.ticks.y=element_blank(),axis.title.y=element_blank()) + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  coord_cartesian(xlim=c(0.5,348.5),ylim=c(0,length(unique(allPathG$Gene))),expand=F,clip="off") + 
  scale_y_continuous(breaks=seq(0.5,length(unique(allPathG$Gene))-0.5,by=1),labels=c(rev(allGene),'Risk Group')) +
  scale_fill_manual(values=fillLeg) + 
  theme(legend.text=element_text(size=15))


grid.newpage()
gt <- ggplotGrob(p)
fig_size <- c(19.5,6) # inches
margin <- unit(2, "line")
gt$vp <- viewport(width = unit(fig_size[1], "in") - margin, height=unit(fig_size[2],"in")- margin)
grid.draw(gt)