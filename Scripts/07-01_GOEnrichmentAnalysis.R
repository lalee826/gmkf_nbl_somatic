##################################################################################################################
### Description: Visualizing enrichment of GO Terms in pathogenic genes
###    

library(parallel)
library(data.table)
library(dplyr)
library(ggplot2)
library(grid)

workdir = '/rocker-build/gmkf_nbl_somatic/'

### read in gene ontology enrichment results
goRes = read.table(paste0(workdir,'Data/GO_analysis.txt'),sep='\t',
                   header=FALSE,stringsAsFactors=FALSE,comment.char='',quote='') %>% 
                  dplyr::select(c('V1','V4','V5','V6','V7','V8'))
colnames(goRes) = c('GO_Process','Expected','Over_Under','Fold_Enrichment','Raw_pVal','FDR_pVAl')
goRes$GO_Process = sapply(goRes$GO_Process,FUN=function(x){trimws(strsplit(x,split='[(]')[[1]][1],which='right')},
                          USE.NAMES=FALSE)


### keep only most specific terms
topPaths = readRDS(paste0(workdir,'Data/topPaths.RDS'))
goRes = goRes %>% filter(GO_Process %in% topPaths)
goRes$FDR_pVAl = as.numeric(goRes$FDR_pVAl)
goRes = goRes[order(goRes$FDR_pVAl),]

### keep GO terms with adjusted p-val <0.001 or if histone related
hisTerms = goRes$GO_Process[sapply(goRes$GO_Process,FUN=function(x){'histone' %in% strsplit(x,split=' ')[[1]]},USE.NAMES=FALSE)]
goResKeep = goRes %>% filter(FDR_pVAl <= 0.001 | GO_Process %in% hisTerms)
#remove some terms we don't care about
goResKeep = goResKeep %>% filter(!(GO_Process %in% c('positive regulation of double-strand break repair',
                                                     'response to estrogen',
                                                     'homeostasis of number of cells within a tissue',
                                                     'fibroblast proliferation',
                                                     'regulation of B cell apoptotic process'
)))

#change some to shorter names for better plots
goResKeep[goResKeep$GO_Process=='positive regulation of transforming growth factor beta receptor signaling pathway',
          'GO_Process'] = "positive regulation of TGFBR signaling pathway"
goResKeep[goResKeep$GO_Process=='negative regulation of transcription by RNA polymerase II','GO_Process'] = "negative regulation of transcription"
goResKeep[goResKeep$GO_Process=='negative regulation of endothelial cell apoptotic process','GO_Process'] = "endothelial cell apoptotic process"
#add color code to terms
redTerms = goResKeep$GO_Process[sapply(goResKeep$GO_Process,FUN=function(x){
  'histone' %in% strsplit(x,split=' ')[[1]] | 'transcription' %in% strsplit(x,split=' ')[[1]]},
  USE.NAMES=FALSE)]
goResKeep$xlabcol = sapply(goResKeep$GO_Process,FUN=function(x){
  ifelse(x %in% redTerms,'red','black')},
  USE.NAMES=FALSE)
goResKeep$GO_Process = factor(goResKeep$GO_Process,levels = goResKeep$GO_Process)

### plot GO terms
p = ggplot() + geom_point(data = goResKeep,aes(x=GO_Process,y=FDR_pVAl),size=4) + theme_minimal() + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.text.x=element_text(size=16,angle=90,color=goResKeep$xlabcol,hjust=0.95,vjust=0.3)) + 
  theme(axis.line.x=element_line(size=2)) + 
  scale_y_continuous(trans = 'log10',breaks=c(0.05,0.02,0.01,0.00001,0.000000001),labels=c(0.05,0.02,0.01,0.00001,0.000000001)) + 
  ylab('FDR Adjusted P-Value') + theme(axis.text.y = element_text(size=15),axis.title.y = element_text(size=18,margin=margin(t=0,r=25,b=0,l=0),angle=0)) +
  theme(axis.title.x = element_blank()) + 
  geom_hline(yintercept=0.001,linetype="dashed",color="gray40",size=1,alpha=0.5)
p

grid.newpage()
gt <- ggplotGrob(p)
panels <- gt$layout$t[grep("panel",gt$layout$name)]
fig_size <- c(18,11) # inches
margin <- unit(2, "line")
gt$vp <- viewport(width = unit(fig_size[1], "in") - margin, height=unit(fig_size[2],"in")- margin)
grid.draw(p)