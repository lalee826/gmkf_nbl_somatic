##################################################################################################################
### Description: Visualizing GSEA of KEGG pathways between ranked lists of pathogenic mutations in GMKF vs.
###              representative subset of pan-pediatric cohort from Grobner et. al (PMID: 29489754)
###             - Grobner cohort consists of:
###               - 30 Acute myeloid leukemias
###               - 20 B-cell acute lymphoblastic leukemia, hypodiploid
###               - 10 B-cell acute lymphoblastic leukemia, non-hypodiploid
###               - 11 Embryonal tumor with multilayered rosettes
###               - 15 Ependymoma infratentorial
###               - 15 Ependymoma supratentorial
###               - 16 Hepatoblastoma
###               - 42 Osteosarcoma
###               - 21 Rhabdomyosarcoma
###               - 19 T-cell acute lymphoblastic leukemia
###               - 10 Medulloblastoma Group3
###               - 10 Medulloblastoma Group4
###               - 10 Medulloblastoma SHH
###               - 10 Medulloblastoma WNT
###               - 25 Pilocytic astrocytoma
###               - 15 "Wilms' tumors
###             - Annotation of pathogenic variants in this external cohort was performed manually as detailed in manuscript
###             - GSEA of KEGG pathways between cohorts was performed using WebGestalt (PMID: 38808672)
###

#import grobner mutation counts
gMutpathGenes = read.table(paste0(workdir,'Data/grobner_path_muts.txt'),
                           header=FALSE,stringsAsFactors=FALSE,sep='\t',comment.char='',quote='')
gMutpathGenes = gMutpathGenes %>% rename('gene'='V1','n'='V2')
#import gmkf mutation counts
kMutpathGenes = read.table(paste0(workdir,'Data/gmkf_path_muts.txt'),
                           header=FALSE,stringsAsFactors=FALSE,sep='\t',comment.char='',quote='')
kMutpathGenes = kMutpathGenes %>% rename('gene'='V1','n'='V2')
#import gmkf gsea results
gGS = read.table(paste0(workdir,'Data/grobner_gsea_results.tsv'),
                 header=TRUE,stringsAsFactors=FALSE,sep='\t',comment.char='',quote='')
#import grobner gsea results
kGS = read.table(paste0(workdir,'Data/gmkf_gsea_results.tsv'),
                 header=TRUE,stringsAsFactors=FALSE,sep='\t',comment.char='',quote='')

#select KEGG pathways to plot
kPath = c('Wnt signaling pathway','Pathways in cancer','AMPK signaling pathway',
          'ATP-dependent chromatin remodeling','Transcriptional misregulation in cancer',
          'NF-kappa B signaling pathway','MicroRNAs in cancer')
gGS$Gene.Set[gGS$Description %in% kPath]
kGS$Gene.Set[kGS$Description %in% kPath]

#need total gene mutation count in these pathways
#import gene lists:
keggGenes = setNames(list(scan(paste0(workdir,'Data/kegg_pathway_genes/hsa04310_genes_cut.txt'),character(),quote=""),
                          scan(paste0(workdir,'Data/kegg_pathway_genes/hsa05200_genes_cut.txt'),character(),quote=""),
                          scan(paste0(workdir,'Data/kegg_pathway_genes/hsa04152_genes_cut.txt'),character(),quote=""),
                          scan(paste0(workdir,'Data/kegg_pathway_genes/hsa03082_genes_cut.txt'),character(),quote=""),
                          scan(paste0(workdir,'Data/kegg_pathway_genes/hsa05202_genes_cut.txt'),character(),quote=""),
                          scan(paste0(workdir,'Data/kegg_pathway_genes/hsa04662_genes_cut.txt'),character(),quote=""),
                          scan(paste0(workdir,'Data/kegg_pathway_genes/hsa04064_genes_cut.txt'),character(),quote=""),
                          scan(paste0(workdir,'Data/kegg_pathway_genes/hsa05206_genes_cut.txt'),character(),quote="")
),
nm=c("hsa04310","hsa05200","hsa04152",
     "hsa03082","hsa05202","hsa04662","hsa04064","hsa05206"))


#make df for plotting
#go through each kegg pathway
totalMutG = gMutpathGenes %>% pull(n) %>% sum()
totalMutK = kMutpathGenes %>% pull(n) %>% sum()
pdf = data.frame(x=c(),y=c(),fType=c(),pval=c())
for (k in kPath) {
  
  #get hsa id,pvalue for kegg pathway
  if (k %in% gGS$Description) {
    hsa = gGS$Gene.Set[gGS$Description==k]
    pv = gGS$P.Value[gGS$Description==k]
    keggType = 'Suppressed'
  } else {
    hsa = kGS$Gene.Set[kGS$Description==k]
    keggType = 'Enriched'
    pv = kGS$P.Value[kGS$Description==k]
  }
  
  geneSet = keggGenes[[hsa]]
  
  #get no. mutations in grobner set
  gn = gMutpathGenes %>% filter(gene %in% geneSet) %>% pull(n) %>% sum()
  #get no. mutations in gmkf set
  kn = kMutpathGenes %>% filter(gene %in% geneSet) %>% pull(n) %>% sum()
  
  #get mut count ratio
  if (gn > kn) {
    gRatio = kn/gn} else {
      gRatio = gn/kn
    }
  
  #create data frame
  dfTemp = data.frame(x=gRatio,y=k,fType=keggType,pval=pv)
  
  #append row
  pdf = rbind(pdf,dfTemp)
  
}

#set y as plotting factor
pdf$y = factor(pdf$y,levels=pdf$y)
pdf$fType = factor(pdf$fType,levels=c('Enriched','Suppressed'))
pdf$sizePlot=1-pdf$pval
pdf$y = factor(pdf$y,levels=c("MicroRNAs in cancer","Pathways in cancer","AMPK signaling pathway","Wnt signaling pathway",
                              "ATP-dependent chromatin remodeling","NF-kappa B signaling pathway","Transcriptional misregulation in cancer"
))

p = ggplot() + geom_point(data=pdf,aes(x=x,y=y,col=pval,size=sizePlot)) +
  #scale_color_distiller(palette = "OrRd") + 
  scale_color_gradient(low='orange',high='blue') +
  facet_grid(cols=vars(fType)) + ylab('KEGG Pathway') + xlab('Mutation Ratio') +
  theme(axis.title=element_text(size=25),axis.text.x=element_text(size=14)) + 
  theme(axis.text.y=element_text(size=15)) + 
  scale_size(range=c(4,16)) +
  theme(panel.background=element_rect(fill='gray100')) +
  theme(panel.grid.major=element_line(color='gray95'),panel.grid.minor=element_blank()) +
  theme(strip.text=element_text(size=25)) + 
  theme(strip.background=element_rect(fill='azure2',color='gray50')) +
  theme(panel.border=element_rect(color='gray50',fill=NA))
p

grid.newpage()
gt <- ggplotGrob(p)
panels <- gt$layout$t[grep("panel",gt$layout$name)]
fig_size <- c(18,11) # inches
margin <- unit(2, "line")
gt$vp <- viewport(width = unit(fig_size[1], "in") - margin, height=unit(fig_size[2],"in")- margin)
grid.draw(p)
