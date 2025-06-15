##################################################################################################################
### Description: Calculate tumor mutation burden in coding regions only using a hg38 transcripts bed file from 
###              USCS genome browser to determine coding bases

library(parallel)
library(data.table)
library(dplyr)
workdir = '/rocker-build/gmkf_nbl_somatic/'

## import coding mutations
# The following samples were found to have likely tumor/normal mismatches or low NGS checkmate correlation in a separate analysis
excludeSamp = readRDS(c(paste0(workdir,'Data/excluded_samples.RDS')))

### load all coding mut data ###
codMuts = read.table(paste0(workdir,'Data/coding_mutations_annotated.tsv'),sep='\t',header=TRUE,row.names=NULL,stringsAsFactors=FALSE)
#remove excluded samples
codMuts <- codMuts %>% filter(!(Tumor_Sample_Barcode %in% excludeSamp))

#remove variants for which prediction by classifier was ambiguous and num callers == 2
codMuts <- codMuts %>% filter(!(numCallers == 2 & prediction %in% c('ambiguous','fail')))

## import tumor clinical data table
clinData = read.csv(paste0(workdir,'Data/tumor_clinical_data.tsv'),sep='\t',header=TRUE,row.names=NULL,stringsAsFactors=FALSE)

##get total bases from ucsc transcripts file
#transcripts bed file
bed = read.table(paste0(workdir,'References/hg38_ucsc_transcripts.bed'),sep='\t',stringsAsFactors=FALSE,header=TRUE,row.names=NULL)
#list of coding genes
genes = readRDS(paste0(workdir,'References/hg38_coding_genes.RDS'))
#filter bed file for coding genes
bed = bed %>% filter(Gene %in% genes)
bed$size = bed$End - bed$Start + 1
#filter out large genes that include TITAN and MUC genes
bed = bed %>% filter(size <= 10e4)

##calculate total size of coding genome used in our data
#ensure we don't double count overlapping bases, go through each chromosome
chrs = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
         'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20',
         'chr21','chr22','chrX','chrY')

#function to count exonic bases in a chromosomes
calc_cod_size = function(chr) {
  
  #print progress if desired
  print(chr)
  
  #filter for chromosome
  bedc = bed %>% filter(Chromosome == chr)
  #calculate in chunks
  s = round(seq(0,nrow(bedc),length.out=10))
  #create vector for all bases in coding genome
  baseInTotal = c()
  
  #iterate through chunks
  for (chunk in c(1:9)) {
    #vector of bases for chunk
    baseIn = c()
    for (i in c((s[chunk]+1):s[chunk+1])) {
      #print(i)
      baseIn = c(baseIn,c(bedc[i,'Start']:bedc[i,'End']))
    }
    #remove duplicates
    baseIn = unique(baseIn)
    baseInTotal = c(baseInTotal,baseIn)
  }
  #remove duplicates
  baseInTotal = unique(baseInTotal)
  #length of vector is size of coding genome
  return(length(baseIn))
}

#count total bases by chromosome
tb = lapply(as.list(chrs),calc_cod_size)
#multi-core operation if desired
#tb = mclapply(as.list(chrs),calc_cod_size,mc.cores = detectCores())

#add total bases and conver to Megabase
tb = sum(unlist(tb)) / 10e5 ##44.8621 MB

###calculate TMB using total bases
#Use vaf cutoff
vcut = 0.125

#sample names
sampleNames = clinData$Tumor_Sample_Barcode
#get mutations per sample
mutBySample = sapply(sampleNames,FUN=function(x){
  codMuts %>% filter(Tumor_Sample_Barcode == x & vaf >= vcut) %>% 
    tally() %>% pull(n) %>% as.integer()
},USE.NAMES=TRUE) 
#calculate TMB
tmbBySample = mutBySample / tb

#calculate median TMB in high-risk samples
hrSamp = clinData %>% filter(Risk == 'high risk') %>% pull(Tumor_Sample_Barcode)
hrMedianTMB = median(tmbBySample[hrSamp])
#calculate median TMB in intermediate-risk samples
intSamp = clinData %>% filter(Risk == 'intermediate risk') %>% pull(Tumor_Sample_Barcode)
intrMedianTMB = median(tmbBySample[intSamp])
#calculate median TMB in intermediate-risk samples
lowSamp = clinData %>% filter(Risk == 'low risk') %>% pull(Tumor_Sample_Barcode)
lrMedianTMB = median(tmbBySample[lowSamp])
