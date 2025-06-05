# This is an example performing with five samples 

library(maftools)
library(dplyr)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(deconstructSigs)
library(VariantAnnotation)

workdir = '/rocker-build/gmkf_nbl_somatic/'


## pull mutation stats from MAF files. 
samples = list.files(paste0(workdir,'Raw Data/Mutect MAF/'))
#remove _mutect tag
samples = sapply(samples,FUN=function(x){paste(strsplit(x,split='_')[[1]][1:2],collapse='_')},USE.NAMES=FALSE)

#gather variant calls and metrics needed for signature analysis
ml = list()
for (i in 1:length(samples)) {
  
  sa = samples[i]
  
  m = read.maf(paste0(workdir,'Raw Data/Mutect MAF/',sa,'_mutect.maf'))
  m1 = m@data
  m2 = m@maf.silent
  m1 = m1[,c("Chromosome","Start_Position","End_Position","Variant_Type","Reference_Allele","Tumor_Seq_Allele1",
             "Tumor_Seq_Allele2","Tumor_Sample_Barcode","Strand","t_depth","t_ref_count","t_alt_count","n_depth","n_ref_count",
             "n_alt_count","Consequence")]
  m1 = m1 %>% mutate(vaf = (t_alt_count)/(t_depth))
  m1 = m1 %>% filter(vaf >= 0.1 & Variant_Type == 'SNP')
  m2 = m2[,c("Chromosome","Start_Position","End_Position","Variant_Type","Reference_Allele","Tumor_Seq_Allele1",
             "Tumor_Seq_Allele2","Tumor_Sample_Barcode","Strand","t_depth","t_ref_count","t_alt_count","n_depth","n_ref_count",
             "n_alt_count","Consequence")]
  m2 = m2 %>% mutate(vaf = (t_alt_count)/(t_depth))
  m2 = m2 %>% filter(vaf >= 0.1 & Variant_Type == 'SNP')
  m = rbind(m1,m2)
  m$caller = rep('m',nrow(m))
  
  #gather lancet variants
  l = read.maf(paste0(workdir,'Raw Data/Lancet MAF/',sa,'_lancet.maf'))
  l1 = l@data
  l2 = l@maf.silent
  l1 = l1[,c("Chromosome","Start_Position","End_Position","Variant_Type","Reference_Allele","Tumor_Seq_Allele1",
             "Tumor_Seq_Allele2","Tumor_Sample_Barcode","Strand","t_depth","t_ref_count","t_alt_count","n_depth","n_ref_count",
             "n_alt_count","Consequence")]
  l1 = l1 %>% mutate(vaf = (t_alt_count)/(t_depth))
  l1 = l1 %>% filter(vaf >= 0.1 & Variant_Type == 'SNP')
  l2 = l2[,c("Chromosome","Start_Position","End_Position","Variant_Type","Reference_Allele","Tumor_Seq_Allele1",
             "Tumor_Seq_Allele2","Tumor_Sample_Barcode","Strand","t_depth","t_ref_count","t_alt_count","n_depth","n_ref_count",
             "n_alt_count","Consequence")]
  l2 = l2 %>% mutate(vaf = (t_alt_count)/(t_depth))
  l2 = l2 %>% filter(vaf >= 0.1 & Variant_Type == 'SNP')
  l = rbind(l1,l2)
  l$caller = rep('l',nrow(l))
  
  #gather strelka variants
  s = read.maf(paste0(workdir,'Raw Data/Strelka MAF/',sa,'_strelka.maf'))
  s1 = s@data
  s2 = s@maf.silent
  s1 = s1[,c("Chromosome","Start_Position","End_Position","Variant_Type","Reference_Allele","Tumor_Seq_Allele1",
             "Tumor_Seq_Allele2","Tumor_Sample_Barcode","Strand","t_depth","t_ref_count","t_alt_count","n_depth","n_ref_count",
             "n_alt_count","Consequence")]
  s1 = s1 %>% mutate(vaf = (t_alt_count)/(t_depth))
  s1 = s1 %>% filter(vaf >= 0.1 & Variant_Type == 'SNP')
  s2 = s2[,c("Chromosome","Start_Position","End_Position","Variant_Type","Reference_Allele","Tumor_Seq_Allele1",
             "Tumor_Seq_Allele2","Tumor_Sample_Barcode","Strand","t_depth","t_ref_count","t_alt_count","n_depth","n_ref_count",
             "n_alt_count","Consequence")]
  s2 = s2 %>% mutate(vaf = (t_alt_count)/(t_depth))
  s2 = s2 %>% filter(vaf >= 0.1 & Variant_Type == 'SNP')
  s = rbind(s1,s2)
  s$caller = rep('s',nrow(s))
  
  #gather vardict variants
  v = read.maf(paste0(workdir,'Raw Data/Vardict MAF/',sa,'_vardict.maf'))
  v1 = v@data
  v2 = v@maf.silent
  v1 = v1[,c("Chromosome","Start_Position","End_Position","Variant_Type","Reference_Allele","Tumor_Seq_Allele1",
             "Tumor_Seq_Allele2","Tumor_Sample_Barcode","Strand","t_depth","t_ref_count","t_alt_count","n_depth","n_ref_count",
             "n_alt_count","Consequence")]
  v1 = v1 %>% mutate(vaf = (t_alt_count)/(t_depth))
  v1 = v1 %>% filter(vaf >= 0.1 & Variant_Type == 'SNP')
  v2 = v2[,c("Chromosome","Start_Position","End_Position","Variant_Type","Reference_Allele","Tumor_Seq_Allele1",
             "Tumor_Seq_Allele2","Tumor_Sample_Barcode","Strand","t_depth","t_ref_count","t_alt_count","n_depth","n_ref_count",
             "n_alt_count","Consequence")]
  v2 = v2 %>% mutate(vaf = (t_alt_count)/(t_depth))
  v2 = v2 %>% filter(vaf >= 0.1 & Variant_Type == 'SNP')
  v = rbind(v1,v2)
  v$caller = rep('v',nrow(v))
  
  #join all callers
  x = rbind(m,s,l,v)
  
  #group by identical calls
  x = x %>% group_by(Chromosome,Start_Position,End_Position,Reference_Allele,Tumor_Seq_Allele1,
                     Tumor_Seq_Allele2,Strand,Tumor_Sample_Barcode) %>% 
    summarise(callers=toString(caller),t_depth=toString(t_depth),t_alt_count=toString(t_alt_count)) %>% 
    as.data.frame()
  #add number of callers
  x$noCallers = sapply(x$callers,FUN=function(x){length(strsplit(x,split=',')[[1]])})
  x$Tumor_Sample_Barcode = rep(sa,nrow(x))
  x = x %>% filter(noCallers > 1)
  
  ml[[i]] = x
}

#### Perform signature analysis using the DeconstructSigs package
# A file of all samples is available to import in the data folder of the repository named SNVlistForSignatureAnalysis.RDS
# ml = readRDS(paste0(workdir,'Data/SNVlistForSignatureAnalysis.RDS'))

#for each variant get average vaf and standardize depth at 100; also add max alt allele detected and filter vars <5
expand_df <- function(x) {
  
  ta100 = c()
  taMaxV = c()
  
  for (i in 1:nrow(x)) {
    
    print(i)
    
    #get allele counts, avg vaf of all callers, and max alt alleles detected
    tdi = sapply(strsplit(x[i,'t_depth'],split=',')[[1]],FUN=function(x){as.numeric(gsub(' ','',x))},USE.NAMES=FALSE)
    tai = sapply(strsplit(x[i,'t_alt_count'],split=',')[[1]],FUN=function(x){as.numeric(gsub(' ','',x))},USE.NAMES=FALSE)
    ta100 = c(ta100,round(mean(tai/tdi),2)*100)
    taMaxV = c(taMaxV,max(tai))
  }
  
  x$tumor_alt_normalized = ta100
  x$tumor_depth_normalized = rep(100,nrow(x))
  x$tumor_ref_normalized = x$tumor_depth_normalized - x$tumor_alt_normalized
  x$max_alt_allele = taMaxV
  x = x %>% filter(max_alt_allele >= 5)
  
  return(x)
  
}

mll = lapply(ml,FUN=expand_df)
#for multi-core
#ml = mclapply(ml,FUN=expand_df,mc.cores=detectCores())

#change field names, keep snps only and substitutions involving A,C,G,T
alldat = data.table::rbindlist(mll,use.names=TRUE,fill=TRUE)
names(alldat)[which(names(alldat)=='Reference_Allele')] = 'REF'
names(alldat)[which(names(alldat)=='Tumor_Seq_Allele2')] = 'ALT'
alldat = alldat %>% filter(REF %in% c('A','C','G','T')) #3,866,746 --> 3,866,746
alldat = alldat %>% filter(ALT %in% c('A','C','G','T')) #3,866,746 --> 3,866,746
mll <- split(alldat,f=alldat$Tumor_Sample_Barcode,drop=FALSE)

#get genome version
genome <- getBSgenome(BSgenome.Hsapiens.UCSC.hg38,masked=FALSE)
seqlevelsStyle(genome) <- 'NCBI'

#make a vrange object for all samples
makevRange <- function(df,in_keep.extra.columns = TRUE,
                       in_seqinfo = NULL, in_seqnames.field = "Chromosome",
                       in_start.field = "Start_Position", in_end.field = "End_Position", in_PID.field = "Tumor_Sample_Barcode",
                       in_strand.field = "Strand",in_subgroup.field='None',
                       verbose_flag=0) {
  
  out_vr <- NULL
  name_list <- tolower(names(df))
  match_list <- c("ref","alt",tolower(in_seqnames.field),
                  tolower(in_start.field),tolower(in_end.field))
  if(length(which(match_list %in% name_list)) == length(match_list)) {
    my_gr <- GenomicRanges::makeGRangesFromDataFrame(df,
                                                     keep.extra.columns=in_keep.extra.columns,
                                                     seqinfo=in_seqinfo,
                                                     seqnames.field=in_seqnames.field,
                                                     start.field=in_start.field,
                                                     end.field=in_end.field,
                                                     strand.field=in_strand.field)
    
    if(!("+" %in% unique(strand(my_gr))) & !("-" %in% unique(strand(my_gr)))){
      strand(my_gr) <- "+"
    }    
    out_vr <- VRanges(seqnames=seqnames(my_gr),ranges=ranges(my_gr),
                      ref=my_gr$REF,alt=my_gr$ALT,totalDepth = my_gr$tumor_depth_normalized, refDepth = my_gr$tumor_ref_normalized, altDepth = my_gr$tumor_alt_normalized)
    if(tolower(in_PID.field) %in% name_list) {
      column_ind <- 
        min(which(tolower(names(mcols(my_gr)))==tolower(in_PID.field)))
      out_vr$PID <- mcols(my_gr)[,column_ind]
      sampleNames(out_vr) <- mcols(my_gr)[,column_ind]
    } else if("pid" %in% name_list) {
      out_vr$PID <- mcols(my_gr)[,"PID"]
      sampleNames(out_vr) <- mcols(my_gr)[,column_ind]
    } else {
      out_vr$PID <- "dummy_PID"
      sampleNames(out_vr) <- "dummyPID"
    }
    if(tolower(in_subgroup.field) %in% name_list) {
      column_ind <- 
        min(which(tolower(names(mcols(my_gr)))==tolower(in_subgroup.field)))
      mcols(out_vr)[,in_subgroup.field] <- mcols(my_gr)[,column_ind]
    } else if("subgroup" %in% name_list) {
      if(verbose_flag==1){
        cat(paste0("YAPSA:::makeVRangesFromDataFrame::warning:",
                   "in_subgroup.field not a valid column name, ",
                   "but default is valid. Retrieving subgroup information.\n"))
      }
      mcols(out_vr)[,in_subgroup.field] <- my_gr$subgroup
    } else {
      if(verbose_flag==1){cat("YAPSA:::makeVRangesFromDataFrame::warning:",
                              "subgroup information missing. ",
                              "Filling up with dummy entries.\n");}
      mcols(out_vr)[,in_subgroup.field] <- "dummy_subgroup"
    }
    seqlengths(out_vr) <- seqlengths(my_gr)    
  } else {
    if(verbose_flag==1){
      cat("YAPSA:::makeVRangesFromDataFrame::error:mismatch in column names, ",
          "return NULL.\n")
    }
  }
  return(out_vr)
}

#use multi-core
allvr <- mclapply(mll,FUN=makevRange,mc.cores=detectCores())


############################# Signature Decomposition #####################################
#build df from vranges for input into deconstruct sigs
ref_list <- unlist(mclapply(allvr,FUN=function(x){x@ref},mc.cores=detectCores()),use.names=FALSE)
alt_list <- unlist(mclapply(allvr,FUN=function(x){x@alt},mc.cores=detectCores()),use.names=FALSE)
pos_list <- unlist(mclapply(allvr,FUN=function(x){x@ranges@start},mc.cores=detectCores()),use.names=FALSE)
chrom_list1 <- unlist(mclapply(allvr,FUN=function(x){x@seqnames@values},mc.cores=detectCores()),use.names=FALSE)
chrom_list2 <- unlist(mclapply(allvr,FUN=function(x){x@seqnames@lengths},mc.cores=detectCores()),use.names=FALSE)
names_list <- unlist(mclapply(allvr,FUN=function(x){x@sampleNames},mc.cores=detectCores()),use.names=FALSE)
chromosome_list <- unlist(sapply(1:length(chrom_list1),FUN=function(i){
  rep(chrom_list1[i],chrom_list2[i])
}))
mut_df <- data.frame("Sample" = names_list, "chr" = chromosome_list, "pos" = pos_list, "ref" = ref_list, "alt" = alt_list)

#create deconstrutsigs object
dSig <- mut.to.sigs.input(mut_df,bsg = BSgenome.Hsapiens.UCSC.hg38)

#import reference mutational signatures matrix. We use COSMIC SBS V3 which is provided in the references folder.
# signatures.ref Either a data frame or location of signature text file, where rows are signatures,
# columns are trinucleotide contexts
cosmicv3 <- read.table(paste0(workdir,'References/sigProfiler_SBS_signatures_v3.csv'),sep=',',header=TRUE,row.names=NULL,stringsAsFactors=FALSE)
sigMat <- t(read.table(paste0(workdir,'References/COSMIC_v3.2_SBS_GRCh38.txt'),sep='\t',header=TRUE,row.names=NULL,stringsAsFactors=FALSE))
colnames(sigMat) = sigMat['Type',]
sigMat = sigMat[2:length(rownames(sigMat)),]
sigMat = sigMat[,colnames(dSig)]

rn = rownames(sigMat)
sigMat = apply(sigMat,MARGIN=2,as.numeric)
rownames(sigMat) = rn
sigMat = data.frame(sigMat)
colnames(sigMat) = colnames(dSig)

#columns should match between sample data and reference data
dim(sigMat)
dim(dSig)

#create a dSig matrix for each patient then bind into one df
# samples <- clusters$labels
samples <- unique(alldat$Tumor_Sample_Barcode)
cosmicsigs <- do.call(rbind,mclapply(as.list(samples),FUN=function(x){whichSignatures(tumor.ref=dSig,sample.id=x,contexts.needed=TRUE,
                                                                                      signatures.ref=sigMat,signature.cutoff=0.02)$weights},mc.cores=12))
#sum of each row (sample) should be close to 1.00
apply(cosmicsigs,MARGIN=1,sum)

#the cosmicsigs object represents estimated mutational signature composition for each tumor
#downstream uses include mutational signature clustering and subsequent visualization in figure 4 of the manuscript