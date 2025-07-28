#* 1. Detect overlaps with any exons of a bed file for BND
############################################################################################################################################
############################################################################################################################################
detect_overlaps_exons_bnd <- function(df,bed,alldata) {
  #*input df format should be a single row df with :
  #*"case_id","type","Gene","Chromosome","Start","End","intersection_type","Events"           
  #*"Breaks1" (list of all breaks on break end 1),"Breaks2 (list of all breaks on break end 2"
  #*pull in breakpoint data from all data file
  
  #filter bed for gene
  g = df[[1,'Gene']]
  bedG = bed[[g]]
  bedG = bedG %>% dplyr::select(Exon,Start,End) %>% rename(Exon_Start = Start,Exon_End = End)
  
  #set tolerance for overlap with exon
  tol = 2
  
  #list to store overlapping intersections
  intList = list()
  listIX = 1
  
  #get chrom of gene
  c = df$Chromosome
  
  #get all breakpoints for translocation on that chromosome
  ev = strsplit(df$Events,split=(','))[[1]]
  cdf = (alldata %>% filter(event %in% ev) %>% as.data.frame())
  pos1e = c()
  #print(cdf)
  for (i in 1:nrow(cdf)) {
    cRow = cdf[i,c('chrom_be1','pos_be1','chrom_be2','pos_be2')]
    #add breakpoints to vector if chromosome matches
    if (cRow[1,'chrom_be1'] == c) {pos1e = c(pos1e,cRow[1,'pos_be1'])}
    if (cRow[1,'chrom_be2'] == c) {pos1e = c(pos1e,cRow[1,'pos_be2'])}
  }
  
  #ensure breakpoints are unique
  pos1e = unique(pos1e)
  
  #remake df with # rows as breakpoints
  df = df %>% slice(rep(1:n(), each = length(pos1e))) %>% dplyr::select(-Breaks2)
  df$Breaks1 = pos1e
  
  #get breakpoints of exons in gene
  pos1t = bedG$Exon_Start - tol
  pos2t = bedG$Exon_End + tol
  
  #convert coordinates to matrices for quick calculaitions
  m1 = matrix(pos1e,nrow=length(pos1e),ncol=length(pos1t)) #matrix of event starts duplicated into columns
  m3 = matrix(rep(pos1t,length(pos1e)),nrow=length(pos1e),byrow=TRUE) #matrix of exon starts
  m4 = matrix(rep(pos2t,length(pos1e)),nrow=length(pos1e),byrow=TRUE) #matrix of exon stops
  
  m3mm1 = m3 - m1 >= 0
  m4mm1 = m4 - m1 >= 0
  m1mm3 = m1 - m3 >= 0
  
  #print(df)
  
  #find cases where bp intersects exon
  case1 = m4mm1 + m1mm3 == 2
  #find row,col indices of true values
  case1IX = which(case1==TRUE,arr.ind=TRUE)
  if (nrow(case1IX) > 0) {
    for (i in 1:nrow(case1IX)) {
      #pull variant data and append gene exon data
      vIX = case1IX[i,1]
      gIX = case1IX[i,2]
      dfTemp = cbind(df[vIX,],bedG[gIX,])
      dfTemp$intersection_type = 'partial'
      intList[[listIX]] = dfTemp
      listIX = listIX + 1
    }
  }
  
  #print(rbindlist(intList))
  return(rbindlist(intList))
  
}


############################################################################################################################################
############################################################################################################################################
#* 2. Detect overlaps with any exons of a bed file for DEL/DUP/INS calls
############################################################################################################################################
############################################################################################################################################
detect_overlaps_exons <- function(df,bed) {
  #input df format should be a single row df with :
  #"case_id","type","Gene","Chromosome","Start","End","intersection_type","Events"           
  #"Breaks1" (list of all breaks on break end 1),"Breaks2 (list of all breaks on break end 2"  
  
  #filter bed for gene
  g = df[[1,'Gene']]
  bedG = bed[[g]]
  bedG = bedG %>% dplyr::select(Exon,Start,End) %>% rename(Exon_Start = Start,Exon_End = End)
  
  #set tolerance for overlap with exon
  tol = 2
  
  #list to store overlapping intersections
  intList = list()
  listIX = 1
  
  #get breakpoints of SVs and all exons in gene
  pos1e = as.integer(unlist(strsplit(df$Breaks1,split=',')))
  pos2e = as.integer(unlist(strsplit(df$Breaks2,split=',')))
  pos1t = bedG$Exon_Start - tol
  pos2t = bedG$Exon_End + tol
  
  #expand df to number of rows of break coordinates
  if (length(pos1e) > 1 | length(pos2e) > 2) {
    maxL = max(length(pos1e),length(pos2e))#get max length of breakends
    df = df %>% slice(rep(1:n(), each = maxL))
    df$Breaks1 = pos1e
    df$Breaks2 = pos2e
  } else {
    maxL = 1
  }
  
  #convert coordinates to matrices for quick calculaitions
  m1 = matrix(pos1e,nrow=maxL,ncol=length(pos1t)) #matrix of event starts duplicated into columns
  m2 = matrix(pos2e,nrow=maxL,ncol=length(pos1t)) #matrix of event stops duplicated into columns
  m3 = matrix(rep(pos1t,maxL),nrow=maxL,byrow=TRUE) #matrix of exon starts
  m4 = matrix(rep(pos2t,maxL),nrow=maxL,byrow=TRUE) #matrix of exon stops
  
  m4mm2 = m4 - m2 >= 0
  m2mm3 = m2 - m3 >= 0
  m3mm1 = m3 - m1 >= 0
  m2mm4 = m2 - m4 >= 0
  m4mm1 = m4 - m1 >= 0
  m1mm3 = m1 - m3 >= 0
  
  #find cases where event ends in exon
  case1 = m3mm1 + m2mm3 + m4mm2 == 3
  #find row,col indices of true values
  case1IX = which(case1==TRUE,arr.ind=TRUE)
  if (nrow(case1IX) > 0) {
    for (i in 1:nrow(case1IX)) {
      #pull variant data and append gene exon data
      vIX = case1IX[i,1]
      gIX = case1IX[i,2]
      if (vIX > nrow(df)) {
        vIX = 1
      }
      dfTemp = cbind(df[vIX,],bedG[gIX,])
      dfTemp$intersection_type = 'partial'
      intList[[listIX]] = dfTemp
      listIX = listIX + 1
    }
  }
  
  #find cases where event starts in exon
  case2 = m2mm4 + m4mm1 + m1mm3 == 3
  #find row,col indices of true values
  case2IX = which(case2==TRUE,arr.ind=TRUE)
  if (nrow(case2IX) > 0) {
    for (i in 1:nrow(case2IX)) {
      #pull variant data and append gene exon data
      vIX = case2IX[i,1]
      gIX = case2IX[i,2]
      if (vIX > nrow(df)) {
        vIX = 1
      }
      dfTemp = cbind(df[vIX,],bedG[gIX,])
      dfTemp$intersection_type = 'partial'
      intList[[listIX]] = dfTemp
      listIX = listIX + 1
    }
  }
  
  #find cases where event surrounds exon
  case3 = m2mm4 + m3mm1 == 2
  #find row,col indices of true values
  case3IX = which(case3==TRUE,arr.ind=TRUE)
  if (nrow(case3IX) > 0) {
    for (i in 1:nrow(case3IX)) {
      #pull variant data and append gene exon data
      vIX = case3IX[i,1]
      gIX = case3IX[i,2]
      if (vIX > nrow(df)) {
        vIX = 1
      }
      dfTemp = cbind(df[vIX,],bedG[gIX,])
      dfTemp$intersection_type = 'subsumes exon'
      intList[[listIX]] = dfTemp
      listIX = listIX + 1
    }
  }
  
  #find cases where event occurs inside exon
  case4 = m4mm2 + m1mm3 == 2
  #find row,col indices of true values
  case4IX = which(case4==TRUE,arr.ind=TRUE)
  if (nrow(case4IX) > 0) {
    for (i in 1:nrow(case4IX)) {
      #pull variant data and append gene exon data
      vIX = case4IX[i,1]
      gIX = case4IX[i,2]
      if (vIX > nrow(df)) {
        vIX = 1
      }
      dfTemp = cbind(df[vIX,],bedG[gIX,])
      dfTemp$intersection_type = 'inside'
      intList[[listIX]] = dfTemp
      listIX = listIX + 1
    }
  }
  
  print(rbindlist(intList))
  return(rbindlist(intList))
  
}

detect_overlaps_exons_atrx <- function(df,bed) {
  #input df format should be a single row df with :
  #"case_id","type","Gene","Chromosome","Start","End","intersection_type","Events"           
  #"Breaks1" (list of all breaks on break end 1),"Breaks2 (list of all breaks on break end 2"  
  
  #filter bed for gene
  g = df[[1,'Gene']]
  bedG = bed[[g]]
  bedG = bedG %>% dplyr::select(Exon,Start,End) %>% rename(Exon_Start = Start,Exon_End = End)
  
  #set tolerance for overlap with exon
  tol = 2
  
  #list to store overlapping intersections
  intList = list()
  listIX = 1
  
  #set breaks of events to min max
  df$Breaks1 = sapply(df$Breaks1,FUN=function(x){
    b = strsplit(x,split=',')[[1]]
    b = sapply(b,FUN=function(y){sub(' ','',y)},USE.NAMES=FALSE)
    return(min(b))
  },USE.NAMES=FALSE)
  df$Breaks2 = sapply(df$Breaks2,FUN=function(x){
    b = strsplit(x,split=',')[[1]]
    b = sapply(b,FUN=function(y){sub(' ','',y)},USE.NAMES=FALSE)
    return(max(b))
  },USE.NAMES=FALSE)
  
  #get breakpoints of SVs and all exons in gene
  pos1e = as.integer(unlist(strsplit(df$Breaks1,split=',')))
  pos2e = as.integer(unlist(strsplit(df$Breaks2,split=',')))
  pos1t = bedG$Exon_Start - tol
  pos2t = bedG$Exon_End + tol
  
  #expand df to number of rows of break coordinates
  if (length(pos1e) > 1 | length(pos2e) > 2) {
    maxL = max(length(pos1e),length(pos2e))#get max length of breakends
    df = df %>% slice(rep(1:n(), each = maxL))
    df$Breaks1 = pos1e
    df$Breaks2 = pos2e
  } else {
    maxL = 1
  }
  
  #convert coordinates to matrices for quick calculaitions
  m1 = matrix(pos1e,nrow=maxL,ncol=length(pos1t)) #matrix of event starts duplicated into columns
  m2 = matrix(pos2e,nrow=maxL,ncol=length(pos1t)) #matrix of event stops duplicated into columns
  m3 = matrix(rep(pos1t,maxL),nrow=maxL,byrow=TRUE) #matrix of exon starts
  m4 = matrix(rep(pos2t,maxL),nrow=maxL,byrow=TRUE) #matrix of exon stops
  
  m4mm2 = m4 - m2 >= 0
  m2mm3 = m2 - m3 >= 0
  m3mm1 = m3 - m1 >= 0
  m2mm4 = m2 - m4 >= 0
  m4mm1 = m4 - m1 >= 0
  m1mm3 = m1 - m3 >= 0
  
  #find cases where event ends in exon
  case1 = m3mm1 + m2mm3 + m4mm2 == 3
  #find row,col indices of true values
  case1IX = which(case1==TRUE,arr.ind=TRUE)
  if (nrow(case1IX) > 0) {
    for (i in 1:nrow(case1IX)) {
      #pull variant data and append gene exon data
      vIX = case1IX[i,1]
      gIX = case1IX[i,2]
      if (vIX > nrow(df)) {
        vIX = 1
      }
      dfTemp = cbind(df[vIX,],bedG[gIX,])
      dfTemp$intersection_type = 'partial'
      intList[[listIX]] = dfTemp
      listIX = listIX + 1
    }
  }
  
  #find cases where event starts in exon
  case2 = m2mm4 + m4mm1 + m1mm3 == 3
  #find row,col indices of true values
  case2IX = which(case2==TRUE,arr.ind=TRUE)
  if (nrow(case2IX) > 0) {
    for (i in 1:nrow(case2IX)) {
      #pull variant data and append gene exon data
      vIX = case2IX[i,1]
      gIX = case2IX[i,2]
      if (vIX > nrow(df)) {
        vIX = 1
      }
      dfTemp = cbind(df[vIX,],bedG[gIX,])
      dfTemp$intersection_type = 'partial'
      intList[[listIX]] = dfTemp
      listIX = listIX + 1
    }
  }
  
  #find cases where event surrounds exon
  case3 = m2mm4 + m3mm1 == 2
  #find row,col indices of true values
  case3IX = which(case3==TRUE,arr.ind=TRUE)
  if (nrow(case3IX) > 0) {
    for (i in 1:nrow(case3IX)) {
      #pull variant data and append gene exon data
      vIX = case3IX[i,1]
      gIX = case3IX[i,2]
      if (vIX > nrow(df)) {
        vIX = 1
      }
      dfTemp = cbind(df[vIX,],bedG[gIX,])
      dfTemp$intersection_type = 'subsumes exon'
      intList[[listIX]] = dfTemp
      listIX = listIX + 1
    }
  }
  
  #find cases where event occurs inside exon
  case4 = m4mm2 + m1mm3 == 2
  #find row,col indices of true values
  case4IX = which(case4==TRUE,arr.ind=TRUE)
  if (nrow(case4IX) > 0) {
    for (i in 1:nrow(case4IX)) {
      #pull variant data and append gene exon data
      vIX = case4IX[i,1]
      gIX = case4IX[i,2]
      if (vIX > nrow(df)) {
        vIX = 1
      }
      dfTemp = cbind(df[vIX,],bedG[gIX,])
      dfTemp$intersection_type = 'inside'
      intList[[listIX]] = dfTemp
      listIX = listIX + 1
    }
  }
  
  print(rbindlist(intList))
  return(rbindlist(intList))
  
}


############################################################################################################################################
############################################################################################################################################
#* 2. Detect overlaps with any exons of a bed file for INV calls
############################################################################################################################################
############################################################################################################################################
detect_overlaps_exons_inv <- function(df,bed) {
  #input df format should be a single row df with :
  #"case_id","type","Gene","Chromosome","Start","End","intersection_type","Events"           
  #"Breaks1" (list of all breaks on break end 1),"Breaks2 (list of all breaks on break end 2"  
  
  #filter bed for gene
  g = df[[1,'Gene']]
  bedG = bed[[g]]
  bedG = bedG %>% dplyr::select(Exon,Start,End) %>% rename(Exon_Start = Start,Exon_End = End)
  
  #set tolerance for overlap with exon
  tol = 2
  
  #list to store overlapping intersections
  intList = list()
  listIX = 1
  
  #get breakpoints of SVs and all exons in gene
  pos1e = as.integer(unlist(strsplit(df$Breaks1,split=',')))
  pos2e = as.integer(unlist(strsplit(df$Breaks2,split=',')))
  pos1t = bedG$Exon_Start - tol
  pos2t = bedG$Exon_End + tol
  
  #if break1 and 2 lengths do not agree, chop off extra
  if (length(pos1e != length(pos2e))) {
    if (length(pos1e) > length(pos2e)) {
      l = length(pos2e)
      pos1e = pos1e[1:l]
    } else if (length(pos2e) > length(pos1e)) {
      l = length(pos1e)
      pos2e = pos2e[1:l]
    }
  }
  
  #expand df to number of rows of break coordinates
  if (length(pos1e) > 1 | length(pos2e) > 2) {
    maxL = max(length(pos1e),length(pos2e))#get max length of breakends
    df = df %>% slice(rep(1:n(), each = maxL))
    df$Breaks1 = pos1e
    df$Breaks2 = pos2e
  } else {
    maxL = 1
  }
  
  #convert coordinates to matrices for quick calculaitions
  m1 = matrix(pos1e,nrow=maxL,ncol=length(pos1t)) #matrix of event starts duplicated into columns
  m2 = matrix(pos2e,nrow=maxL,ncol=length(pos1t)) #matrix of event stops duplicated into columns
  m3 = matrix(rep(pos1t,maxL),nrow=maxL,byrow=TRUE) #matrix of exon starts
  m4 = matrix(rep(pos2t,maxL),nrow=maxL,byrow=TRUE) #matrix of exon stops
  
  m4mm2 = m4 - m2 >= 0
  m2mm3 = m2 - m3 >= 0
  m3mm1 = m3 - m1 >= 0
  m2mm4 = m2 - m4 >= 0
  m4mm1 = m4 - m1 >= 0
  m1mm3 = m1 - m3 >= 0
  
  #find cases where event ends in exon
  case1 = m3mm1 + m2mm3 + m4mm2 == 3
  #find row,col indices of true values
  case1IX = which(case1==TRUE,arr.ind=TRUE)
  if (nrow(case1IX) > 0) {
    for (i in 1:nrow(case1IX)) {
      #pull variant data and append gene exon data
      vIX = case1IX[i,1]
      gIX = case1IX[i,2]
      if (vIX > nrow(df)) {
        vIX = 1
      }
      dfTemp = cbind(df[vIX,],bedG[gIX,])
      dfTemp$intersection_type = 'partial'
      intList[[listIX]] = dfTemp
      listIX = listIX + 1
    }
  }
  
  #find cases where event starts in exon
  case2 = m2mm4 + m4mm1 + m1mm3 == 3
  #find row,col indices of true values
  case2IX = which(case2==TRUE,arr.ind=TRUE)
  if (nrow(case2IX) > 0) {
    for (i in 1:nrow(case2IX)) {
      #pull variant data and append gene exon data
      vIX = case2IX[i,1]
      gIX = case2IX[i,2]
      if (vIX > nrow(df)) {
        vIX = 1
      }
      dfTemp = cbind(df[vIX,],bedG[gIX,])
      dfTemp$intersection_type = 'partial'
      intList[[listIX]] = dfTemp
      listIX = listIX + 1
    }
  }
  
  #find cases where event surrounds exon
  case3 = m2mm4 + m3mm1 == 2
  #find row,col indices of true values
  case3IX = which(case3==TRUE,arr.ind=TRUE)
  if (nrow(case3IX) > 0) {
    for (i in 1:nrow(case3IX)) {
      #pull variant data and append gene exon data
      vIX = case3IX[i,1]
      gIX = case3IX[i,2]
      if (vIX > nrow(df)) {
        vIX = 1
      }
      dfTemp = cbind(df[vIX,],bedG[gIX,])
      dfTemp$intersection_type = 'subsumes exon'
      intList[[listIX]] = dfTemp
      listIX = listIX + 1
    }
  }
  
  #find cases where event occurs inside exon
  case4 = m4mm2 + m1mm3 == 2
  #find row,col indices of true values
  case4IX = which(case4==TRUE,arr.ind=TRUE)
  if (nrow(case4IX) > 0) {
    for (i in 1:nrow(case4IX)) {
      #pull variant data and append gene exon data
      vIX = case4IX[i,1]
      gIX = case4IX[i,2]
      if (vIX > nrow(df)) {
        vIX = 1
      }
      dfTemp = cbind(df[vIX,],bedG[gIX,])
      dfTemp$intersection_type = 'inside'
      intList[[listIX]] = dfTemp
      listIX = listIX + 1
    }
  }
  
  #print(rbindlist(intList))
  return(rbindlist(intList))
  
}




############################################################################################################################################
############################################################################################################################################
#* 4. Detect overlaps with any genes in a bed file for all SV calls
############################################################################################################################################
############################################################################################################################################
#function that finds overlaps with whole gene transcripts
detect_overlaps <- function(df,bed) {
  #note: both input df and bed should be split into list by chromosome
  
  #perform overlap separately for BND and all else
  dfBND = df %>% filter(type == 'BND') 
  df = df %>% filter(type != 'BND')
  
  #perform non-bnd by chromosome
  df = split(df,f = df$chrom_be1)
  dfBNDbe1 = split(dfBND,f = dfBND$chrom_be1)
  dfBNDbe2 = split(dfBND,f = dfBND$chrom_be2)
  
  #set cutoff for how close breakpoint can be to transcript
  nonBNDtol = 2
  
  #list to store nonBND intersections
  cList = list()
  listIX = 1
  
  print('performing intersections for non-bnd SVs...')
  
  for (c in chrs) {
    
    print(c)
    #c = 'chrY'
    dfc = df[[c]]
    bedc = bed[[c]]
    
    pos1e = dfc$pos_be1
    pos2e = dfc$pos_be2
    pos1t = bedc$Start - nonBNDtol
    pos2t = bedc$End + nonBNDtol
    
    m1 = matrix(pos1e,nrow=length(pos1e),ncol=length(pos1t)) #matrix of event starts duplicated into columns
    m2 = matrix(pos2e,nrow=length(pos2e),ncol=length(pos1t)) #matrix of event stops duplicated into columns
    m3 = matrix(rep(pos1t,length(pos1e)),nrow=length(pos1e),byrow=TRUE) #matrix of transcript starts
    m4 = matrix(rep(pos2t,length(pos1e)),nrow=length(pos1e),byrow=TRUE) #matrix of transcript stops
    
    #perform matrix calculations for quick access
    m4mm2 = m4 - m2 >= 0
    m2mm3 = m2 - m3 >= 0
    m3mm1 = m3 - m1 >= 0
    m2mm4 = m2 - m4 >= 0
    m4mm1 = m4 - m1 >= 0
    m1mm3 = m1 - m3 >= 0
    
    #find cases where event ends in transcript
    case1 = m3mm1 + m2mm3 + m4mm2 == 3
    #find row,col indices of true values
    case1IX = which(case1==TRUE,arr.ind=TRUE)
    if (nrow(case1IX) > 0) {
      for (i in 1:nrow(case1IX)) {
        #pull variant data and append gene transcipt data
        vIX = case1IX[i,1]
        gIX = case1IX[i,2]
        dfTemp = cbind(dfc[vIX,],bedc[gIX,])
        dfTemp$intersection_type = 'partial'
        cList[[listIX]] = dfTemp
        listIX = listIX + 1
      }
    }
    
    #find cases where event starts in transcript
    case2 = m2mm4 + m4mm1 + m1mm3 == 3
    #find row,col indices of true values
    case2IX = which(case2==TRUE,arr.ind=TRUE)
    if (nrow(case2IX) > 0) {
      for (i in 1:nrow(case2IX)) {
        #pull variant data and append gene transcipt data
        vIX = case2IX[i,1]
        gIX = case2IX[i,2]
        dfTemp = cbind(dfc[vIX,],bedc[gIX,])
        dfTemp$intersection_type = 'partial'
        cList[[listIX]] = dfTemp
        listIX = listIX + 1
      }
    }
    
    #find cases where event surrounds transcript
    case3 = m2mm4 + m3mm1 == 2
    #find row,col indices of true values
    case3IX = which(case3==TRUE,arr.ind=TRUE)
    if (nrow(case3IX) > 0) {
      for (i in 1:nrow(case3IX)) {
        #pull variant data and append gene transcipt data
        vIX = case3IX[i,1]
        gIX = case3IX[i,2]
        dfTemp = cbind(dfc[vIX,],bedc[gIX,])
        dfTemp$intersection_type = 'subsumes'
        cList[[listIX]] = dfTemp
        listIX = listIX + 1
      }
    }
    
    #find cases where event occurs inside transcript
    case4 = m4mm2 + m1mm3 == 2
    #find row,col indices of true values
    case4IX = which(case4==TRUE,arr.ind=TRUE)
    if (nrow(case4IX) > 0) {
      for (i in 1:nrow(case4IX)) {
        #pull variant data and append gene transcipt data
        vIX = case4IX[i,1]
        gIX = case4IX[i,2]
        dfTemp = cbind(dfc[vIX,],bedc[gIX,])
        dfTemp$intersection_type = 'inside'
        cList[[listIX]] = dfTemp
        listIX = listIX + 1
      }
    }
  }
  
  #bind list into df for non bnd intersections
  dfNonBND = rbindlist(cList)
  rm(cList)
  
  #list to store BND intersections
  cListBNDBE1 = list()
  cListBNDBE1IX = 1
  cListBNDBE2 = list()
  cListBNDBE2IX = 1
  
  #set tolerance for how close breakpoint can be to transcript start/stops
  BNDtol = 5
  
  print('performing intersections for BND breakend 1...')
  
  #analyze breakend 1 first
  for (c in chrs) {
    
    print(c)
    #c = 'chrY'
    dfc = dfBNDbe1[[c]]
    bedc = bed[[c]]
    
    pos1e = dfc$pos_be1
    pos1t = bedc$Start - BNDtol
    pos2t = bedc$End + BNDtol
    
    if (length(pos1e) == 0) {next}
    m1 = matrix(pos1e,nrow=length(pos1e),ncol=length(pos1t)) #matrix of event starts duplicated into columns
    m3 = matrix(rep(pos1t,length(pos1e)),nrow=length(pos1e),byrow=TRUE) #matrix of transcript starts
    m4 = matrix(rep(pos2t,length(pos1e)),nrow=length(pos1e),byrow=TRUE) #matrix of transcript stops
    
    #perform matrix calculations for quick access
    
    m4mm1 = m4 - m1 >= 0
    m1mm3 = m1 - m3 >= 0
    
    #find cases where breakpoint occurs inside transcript
    case1 = m4mm1 + m1mm3 == 2
    #find row,col indices of true values
    case1IX = which(case1==TRUE,arr.ind=TRUE)
    if (nrow(case1IX) > 0) {
      for (i in 1:nrow(case1IX)) {
        #pull variant data and append gene transcipt data
        vIX = case1IX[i,1]
        gIX = case1IX[i,2]
        dfTemp = cbind(dfc[vIX,],bedc[gIX,])
        dfTemp$intersection_type = 'partial'
        cListBNDBE1[[cListBNDBE1IX]] = dfTemp
        cListBNDBE1IX = cListBNDBE1IX + 1
      }
    }
    
    #find cases where breakpoint is within tolerance limit of transcript start
    #find cases where breakpoint is within tolerance limit of transcript end (add later?)
    
  }
  
  print('performing intersections for BND breakend 2...')
  
  #analyze breakend 2
  for (c in chrs) {
    
    print(c)
    #c = 'chr10'
    dfc = dfBNDbe2[[c]]
    bedc = bed[[c]]
    
    pos1e = dfc$pos_be2
    pos1t = bedc$Start
    pos2t = bedc$End
    
    if (length(pos1e) == 0) {next}
    m1 = matrix(pos1e,nrow=length(pos1e),ncol=length(pos1t)) #matrix of event starts duplicated into columns
    m3 = matrix(rep(pos1t,length(pos1e)),nrow=length(pos1e),byrow=TRUE) #matrix of transcript starts
    m4 = matrix(rep(pos2t,length(pos1e)),nrow=length(pos1e),byrow=TRUE) #matrix of transcript stops
    
    #perform matrix calculations for quick access
    
    m4mm1 = m4 - m1 >= 0
    m1mm3 = m1 - m3 >= 0
    
    #find cases where breakpoint occurs inside transcript
    case1 = m4mm1 + m1mm3 == 2
    #find row,col indices of true values
    case1IX = which(case1==TRUE,arr.ind=TRUE)
    if (nrow(case1IX) > 0) {
      for (i in 1:nrow(case1IX)) {
        #pull variant data and append gene transcipt data
        vIX = case1IX[i,1]
        gIX = case1IX[i,2]
        dfTemp = cbind(dfc[vIX,],bedc[gIX,])
        dfTemp$intersection_type = 'partial'
        cListBNDBE2[[cListBNDBE2IX]] = dfTemp
        cListBNDBE2IX = cListBNDBE2IX + 1
      }
    }
  }
  
  #bind bnd lists into df for bnd intersections
  dfBND = rbindlist(append(cListBNDBE1,cListBNDBE2))
  
  #return dataframe of both nonbnd and bnd intersections
  return(rbind(dfNonBND,dfBND))
  
}
