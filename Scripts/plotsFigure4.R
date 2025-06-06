########################## Visualizations - Figure 4
# Description: Code to generate panels for figure 4. Data files are generated from the following scripts:
#

library(parallel)
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(grid)
library(gtable)
workdir = "/rocker-build/gmkf_nbl_somatic/Data/"

### tumor data file
tmb = read.table(paste0(workdir,'tumor_clinical_data.tsv'),sep='\t',stringsAsFactors = FALSE,header=TRUE)

######################################### MUTATION BURDEN PLOT
tmb$group2 <- factor(tmb$group2,levels=c('Low Risk,Non-Amplified','Low Risk,MYCN-Amplified','Intermediate Risk','High Risk,Non-Amplified',
                                         'High Risk,MYCN-amplified'))
sample_order <- tmb$Tumor_Sample_Barcode[order(tmb$group2,tmb$Total_Mutations_vaf_cutoff,decreasing=TRUE)]
tmb$xcord <- sapply(as.character(tmb$Tumor_Sample_Barcode),FUN=function(x){(1:nrow(tmb))[sample_order==x]})
tmb$Tumor_Sample_Barcode <- factor(tmb$Tumor_Sample_Barcode,levels=sample_order)

#add background panel colors
panel_bg <- data.frame(tmb %>% group_by(group2) %>% tally())
panel_bg$group2 <- factor(panel_bg$group2,levels=c('High Risk,MYCN-amplified','High Risk,Non-Amplified','Intermediate Risk','Low Risk,MYCN-Amplified',
                                                   'Low Risk,Non-Amplified'))
panel_bg <- panel_bg[order(panel_bg$group2),]
panel_bg$start <- c(1,(cumsum(panel_bg$n)[1:(nrow(panel_bg)-1)]+1))
panel_bg$end <- c((panel_bg$start-1)[2:nrow(panel_bg)],nrow(tmb))
panel_bg$sample_start <- sapply(panel_bg$start,FUN=function(x){sample_order[x]})
panel_bg$sample_end <- sapply(panel_bg$end,FUN=function(x){sample_order[x]})
panel_bg$sample_start <- factor(panel_bg$sample_start,levels=sample_order)
panel_bg$sample_end <- factor(panel_bg$sample_end,levels=sample_order)
panel_bg$label_pos <- (panel_bg$start + panel_bg$end) / 2
panel_bg$label_pos[4] <- 225
panel_bg$label_height <- c(0.00005,0.00005,0.00005,0.0001,0.00005)
panel_bg$label_graph <- c('High Risk,\nMYCN-amplified','High Risk,\nNon-Amplified','Intermediate\nRisk','Low Risk,\nMYCN-Amplified',
                          'Low Risk,\nNon-Amplified')
panel_bg$label_graph <- factor(panel_bg$label_graph,levels=c('High Risk,\nMYCN-amplified','High Risk,\nNon-Amplified','Intermediate\nRisk','Low Risk,\nMYCN-Amplified',
                                                             'Low Risk,\nNon-Amplified'))
group_meds <- data.frame(tmb %>% group_by(group2) %>% summarise(median=median(TMB.Mut.mb.recalc)))
panel_bg$median <- sapply(panel_bg$group2,FUN=function(x){group_meds[group_meds$group2==x,'median']})

mutCoordsHM <- tmb[,c('Tumor_Sample_Barcode','Has.explanatory.mutations','TMB.Mut.mb.recalc')]
mutCoordsHM$Has.explanatory.mutations <- sapply(as.character(mutCoordsHM$Has.explanatory.mutations),FUN=function(x){ifelse(x=='No','',x)})
mutCoordsHM$xcord <- sapply(as.character(mutCoordsHM$Tumor_Sample_Barcode),FUN=function(x){(1:nrow(tmb))[sample_order==x]})
seg_df <- data.frame(xs=c(225),xe=c(210),ys=c(0.025),ye=c(0.2674863))



###### add pathways to pathogenic mutations
genePathways <- setNames(c('Transcriptional Regulator','Cell Cycle/Apoptosis','Cell Cycle/Apoptosis','Chromatin/epigenome maintence','Chromatin/epigenome maintence','Transcriptional Regulator',
                           'Hypermutant gene','Chromatin/epigenome maintence','Major cell signaling','Major cell signaling','Hypermutant gene','Hypermutant gene','Transcriptional Regulator',
                           'Major cell signaling','Unknown Function','Tumor Suppressor','Major cell signaling','Major cell signaling','Chromatin/epigenome maintence','Major cell signaling','Chromatin/epigenome maintence',
                           'Cell Metabolism','Major cell signaling','Major cell signaling','Major cell signaling','Transcriptional Regulator','ALK','Chromatin/epigenome maintence',
                           'Chromatin/epigenome maintence','Cell Cycle/Apoptosis','Major cell signaling','Major cell signaling','Major cell signaling','Major cell signaling','ALK','ALK','ALK',
                           'ALK','ALK','ALK','ALK','ALK','ALK','ALK','ALK','ALK','ALK','ALK','ALK','ALK','ALK'),
                         nm = c('BS_PKHXM1NB','BS_XH6P7BY8','BS_MA4DE9P2','BS_ZYXTC212','BS_B80Z459E','BS_GAA781GQ','BS_8ZFK0D8Z','BS_C2FYP4KW','BS_H7T26DMD','BS_3KPKZVGX','BS_0TER803R',
                                'BS_YV3MJ2BD','BS_EXDY4JFF','BS_YPKSF86W','BS_VCXG53E7','BS_Y4V9JXG4','BS_F92J3KBA','BS_KQHSSRW3','BS_6KH6CXMH','BS_MCCJ3EMF','BS_9PEJSEZG','BS_87J6X9KJ',
                                'BS_7MXT22QP','BS_WZXHZN97','BS_HQYHVPC9','BS_CH0ZV29W','BS_FJP95JXF','BS_76SQJSHP','BS_QCP0VAZ7','BS_YTQ6Y7FQ','BS_TWN989ZD','BS_VVWYPEHM','BS_4D9SK48F',
                                'BS_E9YW2KHS','BS_NVJJJ7MT','BS_CQC29CK5','BS_MYVQ8YDT','BS_H0MVFRCE','BS_DGV11G69','BS_QZAV7XVF','BS_2R0QBATZ','BS_31MQXA1J','BS_CGGJZ5Q1','BS_BVAVPKKN',
                                'BS_X2D8S1HY','BS_SYQG48V0','BS_CPGZYTHM','BS_MMCZHZE6','BS_1G3QX306','BS_B5MDTA3Q')
)
tmb$pathGenePathway = sapply(as.character(tmb$Tumor_Sample_Barcode),FUN=function(x){genePathways[x]},USE.NAMES=FALSE)
tmb$pathGenePathway = sapply(tmb$pathGenePathway,FUN=function(x){ifelse(is.na(x),'None',x)})
tmb$pathGenePathway <- factor(tmb$pathGenePathway,levels = c("ALK","Transcriptional Regulator","Cell Cycle/Apoptosis","Chromatin/epigenome maintence","Hypermutant gene",
                                                             "Major cell signaling","Unknown Function","Cell Metabolism",
                                                             "None"))
#scale size based on pathogenic mutation
tmb$size = sapply(as.character(tmb$pathGenePathway),FUN=function(x) {ifelse(x=='None','small','big')},USE.NAMES=FALSE)

bg_fills <- brewer.pal(8,'Pastel2')[c(4,6,3,7,5)]
bg_fills[4] <- 'springgreen4'
ptFills <- brewer.pal(9,'Set1')
allFills = c(bg_fills,ptFills)
allFills[length(allFills)] = '#E1E1E1'
mutpermbplot <- ggplot() + #geom_point(data=tmb,alpha=0.5,aes(y=TMB.Mut.mb.,x=xcord,color=pathGenePathway),size=4) + theme_classic() + 
  geom_point(data=tmb,alpha=0.75,aes(y=TMB.Mut.mb.recalc,x=xcord,fill=pathGenePathway,size=size),col='gray10',shape=21) + theme_classic() + 
  scale_size_manual(values=c(6,4),guide='none') +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank()) + 
  theme(axis.text.y=element_text(size=20),axis.title.y=element_text(size=18,margin=margin(t=0,r=25,b=0,l=0),angle=0)) +
  geom_hline(yintercept=median(tmb$TMB.Mut.mb.recalc),linetype="dashed",color="gray40",size=1,alpha=0.5) + 
  geom_hline(yintercept=2,linetype="dashed",color="gray40",size=1,alpha=0.5) +
  geom_hline(yintercept=10,linetype="dashed",color="gray40",size=1,alpha=0.5) +
  annotate(geom="text",x=300,y=0.22,label="Median mutations/MB",color="gray30",size=5) +
  annotate(geom="text",x=300,y=2.45,label="Pediatric High",color="gray30",size=5) +
  annotate(geom="text",x=300,y=13,label="Hypermutated",color="gray30",size=5) +
  ylab('Mutations per MB') +
  theme(axis.line.x=element_blank()) + 
  geom_rect(data=panel_bg,aes(fill=group2),xmin=panel_bg$start-0.5,xmax=panel_bg$end+0.5,ymin=-Inf,ymax=+Inf,alpha=0.25) + 
  scale_fill_manual(values=allFills,breaks = c('High Risk,MYCN-amplified','High Risk,Non-Amplified','Intermediate Risk','Low Risk,MYCN-Amplified',
                                               'Low Risk,Non-Amplified',"ALK","Transcriptional Regulator","Cell Cycle/Apoptosis","Chromatin/epigenome maintence","Hypermutant gene",
                                               "Major cell signaling","Unknown Function","Cell Metabolism",
                                               "None")) + ##### COLOR SCALES ####
geom_label(data=panel_bg,aes(fill=group2,x=label_pos,y=0.02,label=label_graph),size=5,colour='gray30',fontface='bold',angle=90,show.legend=FALSE) +
  geom_segment(data=panel_bg,aes(x=start-0.5,y=median,xend=end+0.5,yend=median,color=group2),size=2.5,alpha=0.75) +
  #
  #geom_segment(data=panel_bg,aes(x=start-0.5,y=median,xend=end+0.5,yend=median,color='black'),size=10,alpha=0.75) +
  scale_color_manual(values=allFills,breaks = c('High Risk,MYCN-amplified','High Risk,Non-Amplified','Intermediate Risk','Low Risk,MYCN-Amplified',
                                                'Low Risk,Non-Amplified',"ALK","Transcriptional Regulator","Cell Cycle/Apoptosis","Chromatin/epigenome maintence","Hypermutant gene",
                                                "Major cell signaling","Unknown Function","Cell Metabolism",
                                                "None"),guide='none') +
  #geom_text(data=mutCoordsHM,aes(x=xcord+14,y=(TMB.Mut.mb.)*1.15,label=Has.explanatory.mutations),color="magenta2",size=4.5) + 
  scale_y_continuous(limits=c(0.01,115),trans='log10',expand=c(0,0),breaks=c(0.01,0.1,1.0,10,100),labels=c(0.01,0.1,1.0,10,100)) + 
  scale_x_continuous(limits=c(0,nrow(tmb)+1),expand=c(0,0)) + 
  geom_segment(data=seg_df,aes(x=xs,xend=xe,y=ys,yend=ye),colour='darkgreen') +
  theme(legend.position="none")
#geom_text(data=tmb_muts,aes(x=xshift,y=yshift,label=gPlot),check_overlap=FALSE,size=4,color='dodgerblue4')
mutpermbplot



##########################################################################################################
##########################################################################################################
##########################################################################################################
#* add in mutation signatures/tumor purity/MES-ADRN plot
## read in mutsig cluster data
msigc = read.table(paste0(workdir,'mutsig_cluster_results.tsv'),stringsAsFactors=FALSE,header=TRUE,row.names=NULL)
tmb$mutsigcluster = sapply(tmb$Tumor_Sample_Barcode,FUN=function(x){msigc[msigc$sample_id==x,'cluster']})
#import mes/adrn scores
maDF = read.table(paste0(workdir,'mes_adrn_scores.tsv'),sep='\t',stringsAsFactors=FALSE,header=TRUE,row.names=NULL)
maDF = maDF[,c('SAMPLE','CLASS')]
colnames(maDF) = c('case_id','feature')
#import gsva quantile results
gDF = read.table(paste0(workdir,'gsva_results_quantiles_plot_allgenes.tsv'),sep='\t',stringsAsFactors=FALSE,header=TRUE,row.names=NULL)
#import purity results
purDF = read.table(paste0(workdir,'tumor_purity_results.tsv'),sep='\t',stringsAsFactors=FALSE,header=TRUE,row.names=NULL)
purDF = purDF[purDF$case_id %in% unique(gDF$case_id),]
#join data
pDF = rbind(maDF,purDF)
#remove samples not in tmb table
allSamp = unique(pDF$case_id)
toRemove = allSamp[!(allSamp %in% tmb$case_id)]
pDF <- pDF %>% filter(!(case_id %in% toRemove))
#get samples with no RNA data
swRNA = unique(pDF$case_id)
swoRNA = tmb$case_id[!(tmb$case_id %in% swRNA)]
#add samples with no data
for (s in swoRNA) {
  for (g in c('purity','mes-adrn')) {
    pDF = rbind(pDF,data.frame(case_id=s,feature=paste0(g,'.','NODATA')))
  }
}
#add in cluster data
pDF = rbind(pDF,data.frame(case_id=tmb[,'case_id'],feature=sapply(tmb$mutsigcluster,FUN=function(x){paste0('Cluster ',x+1)})))
#add y-value for plotting
pDF$val = rep(1,nrow(pDF))
#set factor order
fOrder = c("purity.NODATA","purity.50%","purity.55%","purity.60%",
           "purity.65%","purity.70%","purity.75%","purity.80%","purity.85%","purity.90%","purity.95%",
           "purity.100%",
           'Mesenchymal','Slight Mesenchymal','Slight Adrenergic','Adrenergic','mes-adrn.NODATA',
           "Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5","Cluster 6","Cluster 7","Cluster 8")
pDF$feature = factor(pDF$feature,levels=fOrder)
#set fill colors by feature order
pal = c("gray40","white","#E3EEF8","#CFE1F2","#B5D4E9","#93C4DE","#6BAED6","#4A97C9","#2E7EBB","#1664AB","#084A92","#08306B") #colors for purity
pal = c(pal,'#E4BF47','#F3DD93','#EEADFF','#D424FF','gray40') #colors for mes-adrn status
pal = c(pal,brewer.pal(8,'Set1'))
#set xvalue for each sample
xDict = setNames(tmb$xcord,nm=tmb$case_id)
pDF$xcord = sapply(pDF$case_id,FUN=function(x){xDict[[x]]})

mPlot <- ggplot(data=pDF,aes(x=xcord,y=val,fill=feature)) + geom_bar(stat='identity',width=1) + 
  theme_minimal() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank()) + 
  theme(axis.text.y=element_text(size=13),axis.title.y=element_blank()) + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  scale_fill_manual(values=pal) +
  #scale_fill_manual(values=pal,labels = c('Mesenchymal','Slight Mesenchymal','Slight Adrenergic','Adrenergic','No Data',rep('Negative',5),'Zero',rep('Positive',5),'Strong Outlier','No Data')) +
  geom_hline(yintercept=c(0:3)) +
  theme(axis.line=element_blank()) + 
  theme(legend.margin=margin(t=0,r=0,b=200,l=125)) + 
  geom_vline(xintercept=c(0.5,46.5,104.5,209.5,210.5,348.5)) + 
  coord_cartesian(xlim=c(0.5,348.5),ylim=c(0,3),expand=F,clip="off") + 
  scale_y_continuous(breaks=seq(0.5,2.5,by=1),labels=c('MutSig Cluster','MES/ADRN','tumor purity')) + 
  guides(fill=guide_legend(ncol=1,byrow=TRUE))
mPlot_no_legend = mPlot + theme(legend.position='None')

############# Combine tmb and mutsig plots together
grid.newpage()
gb1 <- ggplotGrob(mutpermbplot)
gb2 <- ggplotGrob(mPlot_no_legend)
gt <- gtable:::rbind_gtable(gb1,gb2,size="last")
panels <- gt$layout$t[grep("panel",gt$layout$name)]
gt$heights[panels] <- unit(c(8,3,3),units='null')
fig_size <- c(18,11) # inches
margin <- unit(2, "line")
gt$vp <- viewport(width = unit(fig_size[1], "in") - margin, height=unit(fig_size[2],"in")- margin)
grid.draw(gt)


##########################################################################################################
##########################################################################################################
##########################################################################################################
############################### Mut sig plots
#stacked bar charts of sig compositions
#remove any cluster with 0 across all samples
csigs = read.table(paste0(workdir,'cosmicv3_sig_contribution.tsv'),stringsAsFactors=FALSE,header=TRUE,row.names=1)
colnames(csigs)
x = csigs %>% summarise(across(SBS1:SBS94,sum)) %>% as.data.frame()
sigRemove = names(x)[x == 0]
csigsR = csigs[,sapply(colnames(csigs),FUN=function(x){!(x %in% sigRemove)},USE.NAMES=FALSE)]
csigsR$case_id = rownames(csigsR)
dim(csigsR) #74 signatures

#create empty df to store data
g = data.frame()

#df needs: case_id,sig,sigvalue,cluster
for (i in 1:nrow(csigsR)) {
  #var = plp[i,'plp_var']
  cid = csigsR[i,'case_id']
  #vc = paste(var,cid,sep=':')
  d = csigsR %>% filter(case_id == cid)
  clus = d[1,'cluster']
  d = d[,sapply(colnames(d),FUN=function(x){(!x %in% c('Tumor_Sample_Barcode','case_id','cluster'))},USE.NAMES=FALSE)]
  #create long form
  dg <- d %>% gather(key=CosmicSignature,value=ratio,colnames(d))
  dg$cluster = rep(clus,nrow(dg))
  dg$case_id = rep(cid,nrow(dg))
  g = rbind(g,dg)
}

#### OPTIONAL:: Convert any signature at less than X% to a value of zero
sigCutoff = 0.05
g$ratio = sapply(g$ratio,FUN=function(x){ifelse(x < sigCutoff,0,x)},USE.NAMES=FALSE)

#remove any signatures that have zero values across all samples
x = g %>% group_by(CosmicSignature) %>% summarise(sum=sum(ratio))
sigRemove = x$CosmicSignature[x$sum == 0]
g = g[sapply(g$CosmicSignature,FUN=function(x){!(x %in% sigRemove)},USE.NAMES=FALSE),]

###plot cases by cluster as stacked bar plots of signature compositions
g$CosmicSignature = factor(g$CosmicSignature,levels = unique(g$CosmicSignature))
g$case_id = factor(g$case_id,levels = sort(unique(g$case_id)))
g$cluster = sapply(g$case_id,FUN=function(x){msigc[msigc$sample_id==x,'cluster']})
g$cluster = g$cluster + 1
g$cluster = as.character(g$cluster)
clustDict = setNames(c('Cluster1','Cluster2','Cluster3','Cluster4','Cluster5','Cluster6','Cluster7','Cluster8'),c(1:8))
g$cluster = sapply(g$cluster,FUN=function(x){clustDict[[x]]},USE.NAMES=FALSE)
g$cluster = factor(g$cluster,levels = unique(sort((g$cluster))))

pal = colorRampPalette(c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(11,"Spectral")))
fill = pal(length(unique(g$CosmicSignature)))
mclustcomp <- ggplot(data=g) + geom_bar(aes(x=case_id,y=ratio,fill=CosmicSignature),position='fill',stat='identity') +
  facet_grid(cols=vars(cluster),scales='free') +
  #facet_wrap(~cluster) +
  theme_minimal() + scale_fill_manual(values=fill) + theme(axis.text.x=element_text(size=10,angle=90,vjust=0.5,hjust=1.1)) + 
  theme(panel.grid.major=element_blank()) + theme(panel.grid.minor=element_blank()) + 
  theme(legend.text=element_text(size=20)) + theme(legend.title=element_text(size=20)) + 
  theme(axis.text.y=element_text(size=22)) + theme(axis.title.y=element_text(size=25,angle=360,margin=margin(t=0,r=10,b=0,l=0))) + 
  theme(axis.title.x=element_blank()) + ylim(c(0,1)) + ylab('Ratio') + theme(strip.text.x = element_text(size = 20)) +
  theme(axis.text.x=element_blank())
mclustcomp_nolegend = mclustcomp + theme(legend.position='None')
grid.newpage()
gp <- ggplotGrob(mclustcomp_nolegend)
fig_size <- c(18,10) # inches
margin <- unit(2, "line")
gp$vp <- viewport(width = unit(fig_size[1], "in") - margin, height=unit(fig_size[2],"in")- margin)
grid.draw(gp)

#### make boxplots of mut burden by cluster
tmb$mutsigcluster = tmb$mutsigcluster + 1
tmb$mutsigcluster = factor(tmb$mutsigcluster,levels=1:8)
tmb$TMB.Mut.mb.recalc_noZero = sapply(tmb$TMB.Mut.mb.recalc,FUN=function(x){ifelse(x==0,0.01,x)})
mclustboxplot <- ggplot(data=tmb) + geom_boxplot(aes(x=mutsigcluster,y=TMB.Mut.mb.recalc_noZero)) + theme_bw() + theme(axis.title.x=element_blank()) + 
  facet_grid(cols=vars(mutsigcluster),scales='free') + theme(strip.text.x = element_text(size = 20)) +
  theme(axis.text.x=element_blank()) + ylab('Tumor Mutation Burden')  +
  scale_y_continuous(trans='log10') + theme(axis.text.y=element_text(size=22)) + 
  theme(axis.title.y=element_text(size=25,angle=360,margin=margin(t=0,r=10,b=0,l=0)))
mclustboxplot

#multiplot
grid.newpage()
gb1 <- ggplotGrob(mclustcomp_nolegend)
gb2 <- ggplotGrob(mclustboxplot)
gt <- gtable:::rbind_gtable(gb1,gb2,size="last")
panels <- gt$layout$t[grep("panel",gt$layout$name)]
gt$heights[panels] <- unit(c(8,1),units='null')
fig_size <- c(18,11) # inches
margin <- unit(2, "line")
gt$vp <- viewport(width = unit(fig_size[1], "in") - margin, height=unit(fig_size[2],"in")- margin)
grid.draw(gt)

#### pie chart by cluster
pal = colorRampPalette(c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(11,"Spectral")))
fill = pal(length(unique(g$CosmicSignature)))
#filter for cluster and summarize contributions
clus = 1
clus = paste0('Cluster',as.character(clus))
g$cluster = as.character(g$cluster)
gc = g %>% filter(cluster == clus)
pg = gc %>% group_by(CosmicSignature) %>% summarise(total = sum(ratio)) %>% as.data.frame()
#sort by total contribution
pg = pg[sort(pg$total,index.return=TRUE)$ix,]
pieLabels = c('SBS 89 (Unkown Etiology)','SBS 58 (Potential Artifact)','SBS37 (Unkown Etiology)','SBS1 (Endogenous Deamination)','SBS 5 (Unkown Etiology)',
              'SBS44 (Defective MMR)','SBS 8 (Unkown Etiology)','SBS 30 (Defective BER)','SBS 32 (Treatment related)','SBS 57 (Potential Artifact)\n',
              'SBS 16 (Unkown Etiology)\n\n','SBS 51 (Potential Artifact)')
pieLabels = c(pieLabels,rep('',(dim(pg)[1])-length(pieLabels)))
pie(pg$total,rev(pieLabels),main='Total Mutation Signature Contributions in Cluster 1',col=fill,cex=1.5,cex.main=2)
####### Cluster 2 - HR enriched
clus = 2
clus = paste0('Cluster',as.character(clus))
gc = g %>% filter(cluster == clus)
pg = gc %>% group_by(CosmicSignature) %>% summarise(total = sum(ratio)) %>% as.data.frame()
#sort by total contribution
pg = pg[sort(pg$total,index.return=TRUE)$ix,]
pieLabels = c('SBS18 (NB-enriched)','SBS40 (Unkown Etiology)','SBS 8 (Unkown Etiology)','SBS39 (Unkown Etiology)','SBS58 (Potential Artifact)',
              'SBS89 (Unkown Etiology)','SBS1 (Endogenous Deamination)','SBS37 (Unkown Etiology)','SBS5 (Unkown Etiology)','SBS32 (Treatment related)')
pieLabels = c(pieLabels,rep('',(dim(pg)[1])-length(pieLabels)))
pie(pg$total,rev(pieLabels),main='Total Mutation Signature Contributions in Cluster 2',col=fill,cex=1.5,cex.main=2)
####### Cluster 3 - HR enriched
clus = 3
clus = paste0('Cluster',as.character(clus))
gc = g %>% filter(cluster == clus)
pg = gc %>% group_by(CosmicSignature) %>% summarise(total = sum(ratio)) %>% as.data.frame()
#sort by total contribution
pg = pg[sort(pg$total,index.return=TRUE)$ix,]
pieLabels = c('SBS18 (NB-enriched)','SBS5 (Unkown Etiology)','SBS58 (Potential Artifact)','SBS8 (Unkown Etiology)',
              'SBS39 (Unkown Etiology)','SBS1 (Endogenous Deamination)','SBS38 (UV Damage)','SBS40 (Unkown Etiology)',
              '\n\nSBS10d (Defectove POLD1)\n\n','\n\n\nSBS3 (HR Deficiency)\n\n')
pieLabels = c(pieLabels,rep('',(dim(pg)[1])-length(pieLabels)))
pie(pg$total,rev(pieLabels),main='Total Mutation Signature Contributions in Cluster 3',col=fill,cex=1.5,cex.main=2)
####### Cluster 4 - uninteristing
clus = 4
clus = paste0('Cluster',as.character(clus))
gc = g %>% filter(cluster == clus)
pg = gc %>% group_by(CosmicSignature) %>% summarise(total = sum(ratio)) %>% as.data.frame()
#sort by total contribution
pg = pg[sort(pg$total,index.return=TRUE)$ix,]
pieLabels = c('SBS39 (Unkown Etiology)','SBS8 (Unkown Etiology)','SBS58 (Potential Artifact)','SBS37 (Unkown Etiology)','SBS18 (NB-enriched)',
              'SBS5 (Unkown Etiology)','SBS40 (Unkown Etiology)','SBS1 (Endogenous Deamination','SBS25 (Treatment Related)')
pieLabels = c(pieLabels,rep('',(dim(pg)[1])-length(pieLabels)))
pie(pg$total,rev(pieLabels),main='Total Mutation Signature Contributions in Cluster 4',col=fill,cex=1.5,cex.main=2)
####### Cluster 5 - TMB related
clus = 5
clus = paste0('Cluster',as.character(clus))
gc = g %>% filter(cluster == clus)
pg = gc %>% group_by(CosmicSignature) %>% summarise(total = sum(ratio)) %>% as.data.frame()
#sort by total contribution
pg = pg[sort(pg$total,index.return=TRUE)$ix,]
pieLabels = c('SBS5 (Unkown Etiology)','SBS1 (Endogenous Deamination)','SBS37  (Unkown Etiology)','SBS32 (Treatment Related)',
              'SBS8 (Unkown Etiology)','SBS58 (Potential Artifact)','SBS54 (Potential Artifact)',
              'SBS18 (NB-enriched)','SBS87 (Treatment Related)','\nSBS92 (BER Deficiency)','\n\nSBS89 (Unkown Etiology)',
              '\n\n\nSBS46 (MMR Deficiency)','\n\n\n\nSBS93 (Unkown Etiology)')
pieLabels = c(pieLabels,rep('',(dim(pg)[1])-length(pieLabels)))
pie(pg$total,rev(pieLabels),main='Total Mutation Signature Contributions in Cluster 5',col=fill,cex=1.5,cex.main=2)
####### Cluster 6
clus = 6
clus = paste0('Cluster',as.character(clus))
gc = g %>% filter(cluster == clus)
pg = gc %>% group_by(CosmicSignature) %>% summarise(total = sum(ratio)) %>% as.data.frame()
#sort by total contribution
pg = pg[sort(pg$total,index.return=TRUE)$ix,]
pieLabels = c('SBS5 (Unkown Etiology)','SBS8 (Unkown Etiology)','SBS89 (Unkown Etiology)','SBS18 (NB-enriched)',
              'SBS1 (Endogenous Deamination)','SBS37  (Unkown Etiology)','SBS 51 (Potential Artifact)',
              'SBS8 (Unkown Etiology)','SBS39 (Unkown Etiology)')
pieLabels = c(pieLabels,rep('',(dim(pg)[1])-length(pieLabels)))
pie(pg$total,rev(pieLabels),main='Total Mutation Signature Contributions in Cluster 6',col=fill,cex=1.5,cex.main=2)
####### Cluster 7
clus = 7
clus = paste0('Cluster',as.character(clus))
gc = g %>% filter(cluster == clus)
pg = gc %>% group_by(CosmicSignature) %>% summarise(total = sum(ratio)) %>% as.data.frame()
#sort by total contribution
pg = pg[sort(pg$total,index.return=TRUE)$ix,]
pieLabels = c('SBS 58 (Potential Artifact)','SBS1 (Endogenous Deamination)','SBS 51 (Potential Artifact)',
              'SBS18 (NB-enriched)','SBS5 (Unkown Etiology)','SBS89 (Unkown Etiology)','SBS44 (Defective MMR)')
pieLabels = c(pieLabels,rep('',(dim(pg)[1])-length(pieLabels)))
pie(pg$total,rev(pieLabels),main='Total Mutation Signature Contributions in Cluster 7',col=fill,cex=1.5,cex.main=2)
####### Cluster 8
clus = 8
clus = paste0('Cluster',as.character(clus))
gc = g %>% filter(cluster == clus)
pg = gc %>% group_by(CosmicSignature) %>% summarise(total = sum(ratio)) %>% as.data.frame()
#sort by total contribution
pg = pg[sort(pg$total,index.return=TRUE)$ix,]
pieLabels = c('SBS40 (Unkown Etiology)','SBS18 (NB-enriched)','SBS89 (Unkown Etiology)','SBS 58 (Potential Artifact)',
              'SBS8 (Unkown Etiology)','SBS1 (Endogenous Deamination)','SBS39 (Unkown Etiology)','SBS37 (Unkown Etiology)')
pieLabels = c(pieLabels,rep('',(dim(pg)[1])-length(pieLabels)))
pie(pg$total,rev(pieLabels),main='Total Mutation Signature Contributions in Cluster 8',col=fill,cex=1.5,cex.main=2)