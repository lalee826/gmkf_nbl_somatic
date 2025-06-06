##################################################################################################################
### Description: Calculate tumor mutation burden in coding regions only using a hg38 transcripts bed file from 
###              USCS genome browser to determine coding bases

library(parallel)
library(data.table)
library(dplyr)
workdir = '/rocker-build/gmkf_nbl_somatic/'

## import coding mutations
# The following samples were found to have likely tumor/normal mismatches or low NGS checkmate correlation in a separate analysis
excludeSamp = c('BS_65F8T6N7','BS_89A4ZCQ1','BS_DHESTKM8','BS_GHWPE3CX','BS_17JDFJ9P','BS_PRMXCJDJ','BS_J7C45G0Y','BS_MTAMZ7KZ')

### load all coding mut data ###
codMuts = read.table(paste0(workdir,'Data/coding_mutations_annotated.tsv'),sep='\t',header=TRUE,row.names=NULL,stringsAsFactors=FALSE)
#remove excluded samples
codMuts <- codMuts %>% filter(!(Tumor_Sample_Barcode %in% excludeSamp))

#remove variants for which prediction by classifier was ambiguous and num callers == 2
codMuts <- codMuts %>% filter(!(numCallers == 2 & prediction %in% c('ambiguous','fail')))

## import tumor clinical data table
tmb = read.csv(paste0(workdir,'Data/tumor_clinical_data.tsv'),sep='\t',header=TRUE,row.names=NULL,stringsAsFactors=FALSE)

##get total bases from ucsc transcripts file
bed = read.table(paste0(workdir,'References/hg38_ucsc_exons.bed'),sep='\t',stringsAsFactors=FALSE,header=TRUE,row.names=NULL) %>% 
  dplyr::select(Gene,Chromosome,Start,End) 
genes = read.table(paste0(workdir,'References/hg38_coding_genes.tsv'),sep='\t',stringsAsFactors=FALSE,header=TRUE) %>% pull(Gene)
#filter bed file for coding genes
dim(bed) #280,195
bed = bed %>% filter(Gene %in% genes)
dim(bed) #245,794
#add gene size
bed$size = bed$End - bed$Start + 1
#remove genes with cds smaller than 300 bp or larger than 15k bp
removeGenes = bed %>% group_by(Gene) %>% summarize(cds=sum(size)) %>% as.data.frame() %>% filter(cds < 300 | cds > 15000) %>% pull(Gene)
bed = bed %>% filter(!(Gene %in% removeGenes))
dim(bed) #232,108

##calculate total size of coding genome used in our data
#ensure we don't double count overlapping bases, go through each chromosome
chrs = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
         'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20',
         'chr21','chr22','chrX','chrY')


#iterate through chromosomes
calcOverlap <- function(chr) {

  print(chr)
  #track total overlapping bases to deleted from total
  tbDel = 0
  
  #filter for chromosome
  bedc = bed %>% filter(Chromosome == chr) %>% dplyr::arrange(Start)
  n = nrow(bedc)
  
  #check overlap of any genes by these cases:
  '
  1. gene1 start inside gene 2
  2. gene2 inside gene 1
  This will capture all overlaps as the bed file is ordered by start position
  '
  
  #iterate through exons by row until 2nd to last row
  for (i in c(1:(n-1))) {

    #indices of all rows from one after current to end
    ixVec = c((i+1):n)
    
    #get start/end of exon being checked
    g1s = bedc[i,'Start']
    g1e = bedc[i,'End']
    
    #only check against exons that start after current exon (bed file is ordered)
    g2s = bedc$Start[(i+1):n]
    g2e = bedc$End[(i+1):n]
    
    #check case 1 (all three operations must be true)
    ixTrue = ixVec[(g2s >= g1s) + (g2s <= g1e) + (g2e > g1e) == 3]
    #if overlaps exist, subtract gene 2 start from gene 1 end
    if (length(ixTrue) > 0) {
      tbDel = tbDel + sum(g1e - bedc$Start[ixTrue])
    }
    
    #check case 2 (all two operations must be true)
    ixTrue = ixVec[(g2s >= g1s) + (g1e >= g2e) == 2]
    #if overlaps exist, take total size of gene 2
    if (length(ixTrue) > 0) {
      tbDel = tbDel + sum(bedc$size[ixTrue])
    }
  }
  return(tbDel)
}

toSubtract = lapply(as.list(chrs),calcOverlap)

#calculate the sum of sizes of all genes in bed file and subtract overlapping bases for total bases
tb = (sum(bed$size) - sum(unlist(toSubtract))) / 10e5 #44.8621 MB


###calculate TMB using total bases
#Use vaf cutoff
vcut = 0.125

#sample names
sampleNames = unique(codMuts$Tumor_Sample_Barcode)
#get mutations per sample
mutBySample = sapply(sampleNames,FUN=function(x){
  codMuts %>% filter(Tumor_Sample_Barcode == x & vaf >= vcut) %>% 
    tally() %>% pull(n) %>% as.integer()
},USE.NAMES=TRUE) 
#calculate TMB
tmbBySample = mutBySample / tb
