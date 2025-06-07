"""
Description: Parse SV calls from Novobreak to create a table of SV events. This script will:
    
    - Join Breakpoints into events with multiple ends
    - Tabulate significant metrics
    - Provide a tabular view of final somatic filtered SV calls made by Novobreak
    - The final events tables can be filtered based on the metrics/counts created in table features
"""

#%% 1. import modules and data
import os,sys
import numpy as np
import pandas as pd
workdir = '/rocker-build/gmkf_nbl_somatic/'

#import clinical data and sample metadata
clinData = pd.read_csv(workdir + 'Data/tumor_clinical_data.tsv',sep='\t',header=0)
metadata = pd.read_table(workdir + '/Data/sample_manifest.csv',sep=',',header=0,index_col=None)

#import normal sample barcodes (identifiers)
pid = [metadata['Kids First Participant ID'].loc[metadata['Kids First Biospecimen ID']==x].reset_index().at[0,'Kids First Participant ID'] for x in clinData['Tumor_Sample_Barcode'].tolist()]
nid = [metadata['Kids First Biospecimen ID'].loc[np.logical_and(metadata['Kids First Participant ID']==x,metadata['sample_type']=='Normal')].reset_index().at[0,'Kids First Biospecimen ID'] for x in pid]
clinData['Normal_Sample_Barcode'] = nid
tumor_normal_conv_table = clinData[['Tumor_Sample_Barcode','Normal_Sample_Barcode']]

#import example file from repo
df = pd.read_table('/rocker-build/gmkf_nbl_somatic/Raw Data/Novobreak VCF/PATBFN.novobreak.clean.vcf',sep='\t',header=0)
chromosomes = ['chr' + str(i) for i in np.arange(1,23)] + ['chrX','chrY']

#%% 5. split the info column
#########################################################
#check columns in vcf
#df.columns.values

###extract data from info and tumor/normal columns
info = df.INFO
quals = df.QUAL
fils = df.FILTER

###create these to make table of sv calls
info_dicts = [] #metrics for variant call in list of dictionaries
is_precise = [] #bool for IMPRECISE metric
is_imprecise = []
is_somatic = []
is_germline = []
is_unklen = []
#is_somatic = [] #bool for SECONDARY metric
#tnv = [] #metrics for tumor and normal in list of dictionaries


for i in np.arange(0,len(info)):
    
    #split info field
    infoSplit = info[i].split(';')
    #split each info field by key,value
    infoSplit = [x.split('=') for x in infoSplit]
    #check for flags: precise,imprecise,somatic,germline,unknown length
    isImprecise = np.any([x[0]=='IMPRECISE' for x in infoSplit])
    isPrecise = np.any([x[0]=='PRECISE' for x in infoSplit])
    isSomatic = np.any([x[0]=='SOMATIC' for x in infoSplit])
    isGermline = np.any([x[0]=='GERMLINE' for x in infoSplit])
    isUnkLen = np.any([x[0]=='UNKNOWN_LEN' for x in infoSplit])
    '''
    #remove FLAGS and check if any other fields are len==1
    if np.any([len(y) != 2 for y in list(np.array(infoSplit)[[x[0] not in ['IMPRECISE','PRECISE','SOMATIC','GERMLINE','UNKOWN_LEN'] for x in infoSplit]])]):
        print(i)
    '''
    #remove any fields that are not of key=value format
    infoSplit = list(np.array(infoSplit,dtype=object)[[len(x)==2 for x in infoSplit]])
    #create dictionary of fields and values
    #save row to lists
    info_dicts.append(dict(tuple([tuple(x) for x in infoSplit])))
    is_imprecise.append(isImprecise)
    is_precise.append(isPrecise)
    is_somatic.append(isSomatic)
    is_germline.append(isGermline)
    is_unklen.append(isUnkLen)

#%% 6. Extract the values to be included in table
######### PRIMARY METRICS
#pull chroms
chroms = df.CHROM.tolist()
#pull POS
pos = df.POS.tolist()
#pull SVTYPE
types = [info_dicts[i]['SVTYPE'] if 'SVTYPE' in info_dicts[i] else 'NaN' for i in np.arange(0,len(info_dicts))]
types = ['BND' if x == 'TRA' else x for x in types]
#give a unique id for each sv
ids = ['nb:'+types[i]+':'+chroms[i]+':'+str(pos[i]) for i in np.arange(0,len(chroms))]
#pull end position
ends = [info_dicts[i]['END'] if 'END' in info_dicts[i] else 'NaN' for i in np.arange(0,len(info_dicts))]
#pull sv length
svlen = [abs(int(info_dicts[i]['SVLEN'])) if 'SVLEN' in info_dicts[i] else 'NaN' for i in np.arange(0,len(info_dicts))]
#pull consensus sequence
cons = [info_dicts[i]['CONSENSUS'] if 'CONSENSUS' in info_dicts[i] else 'NaN' for i in np.arange(0,len(info_dicts))]
#pull mate chromosome
chroms2 = [info_dicts[i]['CHR2'] if 'CHR2' in info_dicts[i] else 'NaN' for i in np.arange(0,len(info_dicts))]
#pull Paired-end signature induced connection type
ct = [info_dicts[i]['CT'] if 'CT' in info_dicts[i] else 'NaN' for i in np.arange(0,len(info_dicts))]
#pull Confidence interval around POS for imprecise variants
cipos = [info_dicts[i]['CIPOS'] if 'CIPOS' in info_dicts[i] else 'NaN' for i in np.arange(0,len(info_dicts))]
#pull Confidence interval around END for imprecise variants
ciend = [info_dicts[i]['CIEND'] if 'CIEND' in info_dicts[i] else 'NaN' for i in np.arange(0,len(info_dicts))]
#create a mate id
mid = ['nb:'+types[i]+':'+chroms2[i]+':'+str(ends[i]) if types[i] == 'BND' else 'NaN' for i in np.arange(0,len(chroms))]

#%% 7. create a data frame to output
svdf = pd.DataFrame({'id':ids,'mateid':mid,'type':types,'chrom':chroms,'chrom_mate':chroms2,'pos':pos,'end':ends,'length':svlen,'quality':quals,'filter':fils, #primary metrics
                     'precise':is_precise,'imprecise':is_imprecise,'somatic':is_somatic,'germline':is_germline,'unknown_length':is_unklen, #sv flags
                     'consensus_sequence':cons,'connection_type':ct,'confidence_interval_pos':cipos,'confidence_interval_end':ciend}) #info tags
#add in rest of novobreak fields
svdf = pd.concat([svdf,df.iloc[:,10:]],axis=1)
#%% 8. combine mate pairs and make events df
#### For BND events join mate pairs, discard those without mates
#split off BND from rest
svdf_bnd = svdf.loc[svdf.type=='BND',:].reset_index(drop=True)
svdf = svdf.loc[svdf.type!='BND',:].reset_index(drop=True)

#%% 9. find bnd svs that have mates
#discard rows whose mate pair is not in data frame
#match by mate chromosome and position
match_tolerance = 1e3 #positions must match within this distance
keep = [] #these are indices to keep (mate pair exists)
totalMatePairs = [] #for each tally the number of mates
whichMatePairs = [] #for each keep track of mate IDS
distanceToPair = [] #for each keep track of distance to pair coordinate
#add new row to keep track if cluster ids need to be combined between mates
newClusterID = [] #for each keep track if cluster ids of row and mate need to be combined

for ix,row in svdf_bnd.iterrows():
    #row = svdf_bnd.loc[346,:]
    #match by chromosomes
    t = svdf_bnd.loc[svdf_bnd.chrom==row['chrom_mate'],:] #mate's chrom matches row
    t = t.loc[t.chrom_mate==row['chrom']] #row's chrom matches mate's chrom
    n = t.shape[0] #total rows matching by chrom
    #get positions and match against end position of current row
    p = t.pos.tolist()
    distanceP = [abs(int(row['end']) - int(x)) for x in p] #distance to each position
    hasMatchP = [x <= match_tolerance for x in distanceP]
    distanceP = np.array(distanceP)[hasMatchP].tolist()
    #get end positions and match against position of current row
    e = t.end.tolist()
    distanceE = [abs(int(row['pos']) - int(x)) for x in e] #distance to each end position
    hasMatchE = [x <= match_tolerance for x in distanceE]
    distanceE = np.array(distanceE)[hasMatchE].tolist()
    
    #both ends must match to be considered a mate pair
    hasMatch = np.logical_and(hasMatchE,hasMatchP)

    #if mate exists save index and indices of matches
    #matchIX = []
    if (np.any(hasMatch)):
        keep.append(ix) #save the index to keep
        t = t.loc[[hasMatch[i] == True for i in np.arange(0,n)]] #keep matching rows
        mateIDS = t.id.tolist()
        #mateIDS = t.id[[hasMatch[i] == True for i in np.arange(0,n)]].tolist()
        whichMatePairs.append(mateIDS)
        #in case mate pairs are from different clusters, save a new field that combines both clusters
        #nClusters = len(set([row['cluster_id']] + svdf_bnd.loc[[x in mateIDS for x in svdf_bnd.id],'cluster_id'].tolist()))
        nClusters = len(set([row['cluster_id']] + t.cluster_id.tolist()))
        #if the mate pairs come from different clusters then join cluster ids
        if nClusters > 1:
            c = list(set([row['cluster_id']] + list(set(t['cluster_id']))))
            c.sort()
            newClusterID.append((':').join(c))
        else: newClusterID.append('NaN')
    else:
        whichMatePairs.append('NaN')
        
#create dictionary to hold cluster id mates   
cluster_mates_dict = dict(zip(svdf_bnd.cluster_id,[[] for i in np.arange(0,len(svdf_bnd.cluster_id))]))
#filter for sv with mates
svdf_bnd = svdf_bnd.loc[keep,:].reset_index(drop=True)
#add list of matepairs
whichMatePairs = np.array(whichMatePairs,dtype=object)[[x != 'NaN' for x in whichMatePairs]].tolist()
svdf_bnd['mates'] = whichMatePairs
#add new cluster id if mates do not match
svdf_bnd['cluster_id'] = [svdf_bnd['cluster_id'][i] if newClusterID[i] == 'NaN' else newClusterID[i] for i in np.arange(0,svdf_bnd.shape[0])]

### if one cluster belongs to multiple new ids, join all ids together
#fill dictionary with all cluster id mates
newClusterID_split = [x.split(':') for x in newClusterID]
for cGroup in newClusterID_split:
    if len(cGroup) == 1: continue
    else:
        c = np.array(cGroup)
        for i in np.arange(0,len(c)):
             cluster_mates_dict[c[i]] = list(set(cluster_mates_dict[c[i]] + c[[k != i for k in np.arange(0,len(c))]].tolist()))

#now add all mate clusters
cid = svdf_bnd['cluster_id'].tolist()
for i in np.arange(0,len(cid)):
    for it in [1,2,3,4,5]: #perform two rounds of joining
        c = cid[i].split(':')
        nc = []
        for x in c:
            nc += cluster_mates_dict[x] + [x]
        nc = list(set(nc))
        nc.sort()
        cid[i] = (':').join(nc)
svdf_bnd['cluster_id'] = cid

#record any new cluster ids
newClusterID = newClusterID + np.array(cid)[[len(x.split(':')) > 1 for x in cid]].tolist()

### we must check that we have all appropriate rows in new cluster ids that were created
for ev in set(newClusterID):
    idChange = [] #ids that need to be changed
    if ev == 'NaN': continue
    else:
        t = svdf_bnd.loc[svdf_bnd.cluster_id==ev,:] #all rows that have combined id
        if t.shape[0] == 0: continue
        d = svdf_bnd.loc[svdf_bnd.cluster_id!=ev,:] #all rows that do not
        #get chroms in event
        chrs = list(set(t.chrom))
        chr1,chr2 = chrs[0],chrs[1]
        #identify any rows that should belong
        dd = d.loc[np.logical_or(np.logical_and(d.chrom==chr1,d.chrom_mate==chr2),np.logical_and(d.chrom==chr2,d.chrom_mate==chr1))]
        #keep any rows that are within 1kb of the breakpoints
        if dd.shape[0] > 0:
            for ix,row in dd.iterrows():
                #idenitfy rows with combined id that have this chromosome orientation
                tt = t.loc[np.logical_and(t.chrom==row['chrom'],t.chrom_mate==row['chrom_mate'])]
                #check if pos matches all in combined id rows
                posMatch = np.all([abs(row['pos']-x) < 1000 for x in tt.pos])
                #chec if end matches all in combined id rows
                endMatch = np.all([abs(int(row['end'])-int(x)) < 1000 for x in tt.end])
                #if both true then event id should be changed
                if posMatch and endMatch: idChange.append(row['id'])
        #change cluster id for all ids that matched
        for id in idChange:
            svdf_bnd.loc[svdf_bnd.id==id,'cluster_id'] = ev

#%% 10. if there are no bnd events, just output the non bnd svs and exit program
if svdf_bnd.shape[0] == 0:
    #~~~~~~~~~~~~~~~~create svdf ready for output
    #create same data frame structure as other samples
    split_cols = ['id','chrom','pos','quality','filter','precise','imprecise','somatic','germline','unknown_length','consensus_sequence','connection_type',
              'confidence_interval_pos','confidence_interval_end','cluster_id','contig_id','contig_size','reads_used_for_assembly','average_coverage',
              'tumor_bkpt1_depth','tumor_bkpt1_sp_reads','tumor_bkpt1_qual','tumor_bkpt1_high_qual_sp_reads','tumor_bkpt1_high_qual_qual',
              'normal_bkpt1_depth','normal_bkpt1_sp_reads','normal_bkpt1_qual','normal_bkpt1_high_qual_sp_reads','normal_bkpt1_high_qual_qual',
              'tumor_bkpt2_depth','tumor_bkpt2_sp_reads','tumor_bkpt2_qual','tumor_bkpt2_high_qual_sp_reads','tumor_bkpt2_high_qual_qual',
              'normal_bkpt2_depth','normal_bkpt2_sp_reads','normal_bkpt2_qual','normal_bkpt2_high_qual_sp_reads','normal_bkpt2_high_qual_qual',
              'tumor_bkpt1_discordant_reads','normal_bkpt1_discordant_reads','tumor_bkpt2_discordant_reads','normal_bkpt2_discordant_reads']
    #modify non bnd columns and names and bind with BND events
    rename_dict = dict(zip(split_cols,[x+'_be1' for x in split_cols]))
    rename_dict['end'] = 'pos_be2' #end coordinate is breakend2 
    rename_dict['chrom_mate'] = 'chrom_be2' #mate chrom is breakend2 chrom
    svdf = svdf.rename(columns=rename_dict)
    #drop mateid
    svdf = svdf.drop(['mateid'],axis=1)
    #add NaN values for be 1 alt and be2 columns
    nanFill = ['NaN']*svdf.shape[0]
    new_cols_dict = dict(zip([q for p in [[x+y for y in ['_be1_alt','_be1_alt2','_be2','_be2_alt','_be2_alt2']] for x in split_cols] for q in p],[nanFill]*(len(split_cols)*5)))
    #don't change pos_be2 or chrom_be2
    del new_cols_dict['pos_be2']
    del new_cols_dict['chrom_be2']
    
    #add event and mates col
    svdf['event'] = [(':').join([svdf.cluster_id_be1[i],svdf.id_be1[i]]) for i in np.arange(0,svdf.shape[0])]
    svdf['mates'] = [svdf.id_be1[i] for i in np.arange(0,svdf.shape[0])]
    
    #add new columns to non bnd structural variants
    svdf = pd.concat([svdf,pd.DataFrame(new_cols_dict)],axis=1)
    #put cols in order
    col_order = ['event','type','mates', #common columns
             'id_be1','id_be1_alt','id_be1_alt2','id_be2','id_be2_alt','id_be2_alt'] + [q for p in [[y+x for y in ['chrom','pos']]
                   for x in ['_be1','_be1_alt','_be1_alt2','_be2','_be2_alt','_be2_alt2']] for q in p] + ['length'] + [q for p in [[y+x for y in ['quality',
              'filter','precise','imprecise','somatic','germline','unknown_length','consensus_sequence','connection_type','confidence_interval_pos',
              'confidence_interval_end','cluster_id','contig_id','contig_size','reads_used_for_assembly','average_coverage','tumor_bkpt1_depth',
              'tumor_bkpt1_sp_reads','tumor_bkpt1_qual','tumor_bkpt1_high_qual_sp_reads','tumor_bkpt1_high_qual_qual','normal_bkpt1_depth',
              'normal_bkpt1_sp_reads','normal_bkpt1_qual','normal_bkpt1_high_qual_sp_reads','normal_bkpt1_high_qual_qual','tumor_bkpt2_depth',
              'tumor_bkpt2_sp_reads','tumor_bkpt2_qual','tumor_bkpt2_high_qual_sp_reads','tumor_bkpt2_high_qual_qual','normal_bkpt2_depth',
              'normal_bkpt2_sp_reads','normal_bkpt2_qual','normal_bkpt2_high_qual_sp_reads','normal_bkpt2_high_qual_qual','tumor_bkpt1_discordant_reads',
              'normal_bkpt1_discordant_reads','tumor_bkpt2_discordant_reads','normal_bkpt2_discordant_reads']] for x in ['_be1','_be1_alt','_be1_alt2',
              '_be2','_be2_alt','_be2_alt2']] for q in p]
    svdf = svdf[col_order] 
    
    ### change data types
    for col in ['chrom_be1','chrom_be1_alt','chrom_be1_alt2','chrom_be2','chrom_be2_alt','chrom_be2_alt2']:
        svdf[col] = pd.Categorical(svdf[col], 
                          categories=chromosomes + ['NaN'],
                          ordered=True)    
    dtype_dict = {'pos_be1':int,'pos_be1_alt':str,'pos_be1_alt2':str,'pos_be2':int,'pos_be2_alt':str,'pos_be2_alt2':str,'length':int}
    '''
    dtype_dict = {'pos_be1':int,'pos_be1_alt':str,'pos_be1_alt2':str,'pos_be2':int,'pos_be2_alt':str,'pos_be2_alt2':str,
                  'precise_be1':bool,'precise_be1_alt':bool,'precise_be1_alt2':bool,'precise_be2':bool,'precise_be2_alt':bool,'precise_be2_alt2':bool,
                  'imprecise_be1':bool,'imprecise_be1_alt':bool,'imprecise_be1_alt2':bool,'imprecise_be2':bool,'imprecise_be2_alt':bool,'imprecise_be2_alt2':bool,
                  'somatic_be1':bool,'somatic_be1_alt':bool,'somatic_be1_alt2':bool,'somatic_be2':bool,'somatic_be2_alt':bool,'somatic_be2_alt2':bool,
                  'germline_be1':bool,'germline_be1_alt':bool,'germline_be1_alt2':bool,'germline_be2':bool,'germline_be2_alt':bool,'germline_be2_alt2':bool,
                  'unknown_length_be1':bool,'unknown_length_be1_alt':bool,'unknown_length_be1_alt2':bool,'unknown_length_be2':bool,'unknown_length_be2_alt':bool,'unknown_length_be2_alt2':bool,
                  'length':int}
    '''
    svdf = svdf.astype(dtype_dict)
    #order by breakend 1 chrom and pos
    svdf = svdf.sort_values(['chrom_be1','pos_be1','chrom_be2','pos_be2']).reset_index(drop=True)
    
    #save output
    svdf.to_csv(outPath,sep='\t',header=True,index=False)
    
    #exit program
    sys.exit("No BND with mates; saved non BND SVs for sample "+sample+". Exiting now...")

#%% 11. group bnd by event
'''
Events are identified by the cluster id
Need to first check for clusters that have more than one event
For multiple events in same cluster ID separate into multiple events
'''
svdf_bnd['event'] = ['NaN']*svdf_bnd.shape[0]

for ix,ev in enumerate(set(svdf_bnd.cluster_id)):
    '''
    ev = 'E129152:E181081:E234783:E30391:E360502:E362720:E372229:E372798:E394994:E429634:E549137:E555055:E572861:E591950:E629337:E717421:E762237:E8894'
    '''
    t = svdf_bnd.loc[svdf_bnd.cluster_id==ev,].reset_index()
    t.end = t.end.astype(int)
    n = t.shape[0]
    #find all breakpoint clusters
    bc = []
    for i in np.arange(0,n):
        #if first row, add cluster
        if i == 0:
            bc.append({'chrom':t.loc[i,'chrom'],'pos':{t.loc[i,'pos']},'ix':[i],'chrom_mate':t.loc[i,'chrom_mate'],'pos_mate':{t.loc[i,'end']}})
        else:
            #check each saved cluster
            newCluster = True
            for j in np.arange(0,len(bc)):
                #first check chrom and mate chrom
                if t.loc[i,'chrom'] == bc[j]['chrom'] and t.loc[i,'chrom_mate'] == bc[j]['chrom_mate']:
                    #check each position
                    if np.all([np.any([abs(t.loc[i,'pos'] - x) <= 1200 for x in bc[j]['pos']]),
                               np.any([abs(t.loc[i,'end'] - x) <= 1200 for x in bc[j]['pos_mate']])]):
                        bc[j]['pos'] = set(list(bc[j]['pos']) + [t.loc[i,'pos']]) #add new position to cluster
                        bc[j]['pos_mate'] = set(list(bc[j]['pos_mate']) + [t.loc[i,'end']]) #add new mate positions to cluster
                        bc[j]['ix'] = bc[j]['ix'] + [i]#add new index to cluster
                        newCluster = False #set new cluster flag to false
                        break
                else: 
                    continue
            if newCluster: bc.append({'chrom':t.loc[i,'chrom'],'pos':{t.loc[i,'pos']},'ix':[i],'chrom_mate':t.loc[i,'chrom_mate'],
                                      'pos_mate':{t.loc[i,'end']}})
    #check if any breakpoint clusters should be joined by checking if chroms match and any position is w/in 1kb of any in other
    bcNew = []
    skipSet = []
    for i in np.arange(0,len(bc)):
        if i in skipSet: continue
        b = bc[i]    
        for j in np.arange(i+1,len(bc)):        
            if j in skipSet: continue
            #check chrom and mate chrom match/check pos and mate pos
            if bc[i]['chrom'] == bc[j]['chrom'] and bc[i]['chrom_mate'] == bc[j]['chrom_mate']: 
                #check each position and mate position
                if np.all([np.any([q for p in [[abs(x-y) <= 1200 for y in bc[i]['pos']] for x in bc[j]['pos']] for q in p]),
                    np.any([q for p in [[abs(x-y) <= 1200 for y in bc[i]['pos_mate']] for x in bc[j]['pos_mate']] for q in p])]):
                    #add cluster j to i
                    b['pos'] = set(list(b['pos']) + list(bc[j]['pos']))
                    b['pos_mate'] = set(list(b['pos_mate']) + list(bc[j]['pos_mate']))
                    b['ix'] = list(set(list(b['ix']) + bc[j]['ix']))
                    skipSet.append(j)
        bcNew.append(b)

    #for each breakpoint cluster, find parter cluster
    tchroms = t.chrom #list of chromosomes
    tpos = t.pos #list of positions
    mateChroms = t.chrom_mate #list of mate chromosomes
    matePos = t.end.astype(int) #list of mate positions
    mateSets = set() #save each set of indices in a event
    for i in np.arange(0,len(bcNew)):
        partners = np.arange(0,n)[[tchroms[j] == bcNew[i]['chrom_mate'] and mateChroms[j] == bcNew[i]['chrom'] and 
                                   np.any([abs(matePos[j] - k) <= 1200 for k in list(bcNew[i]['pos'])]) and
                                   np.any([abs(tpos[j] - k) <= 1200 for k in list(bcNew[i]['pos_mate'])]) for j in np.arange(0,n)]]
        bcNew[i]['mates'] = partners.tolist()
        allIX = (bcNew[i]['ix']+bcNew[i]['mates'])
        allIX.sort()
        mateSets.add(tuple(allIX))
    #give a unique event id for each mate set
    setNo = 1
    for ms in mateSets:
        eventid = ev+'-'+str(setNo)
        t.loc[list(ms),'event'] = eventid
        t.loc[list(ms),'mates'] = ('|').join(t.loc[list(ms),'id'].tolist())
        setNo += 1
    t = t.set_index('index')
    svdf_bnd.loc[t.index.values,'event'] = t['event']
    svdf_bnd.loc[t.index.values,'mates'] = t['mates']

'''
for each event, ensure that there are:
    - > 1 SVs
    - <= 2 chromosomes
    - if two chromosomes: one break location on each chromosome (all breaks w/in 1kb of each other)
    - if one chromosome: two break locations for which each has breaks within 2kb of each other
'''
for ix,ev in enumerate(set(svdf_bnd.event)):
    #ensure more than 1 sv
    t = svdf_bnd.loc[svdf_bnd.event==ev,].reset_index(drop=True)
    numSV = t.shape[0]
    try:
        assert numSV > 1
    except: raise AssertionError('Sample: ' + sample + ';There is only 1 sv in the event ',ev)
    #ensure there are no more than 2 chromosomes in event
    numChrom = set(t.chrom)
    try:
        assert len(numChrom) <= 2
    except: raise AssertionError('Sample: ' + sample + ';There are more than 2 chromosomes in the event ',ev)
    #if there are two chromosomes check pos of each chromosome to check each are <= 1kb of each other
    if len(numChrom) == 2:
        g = t.groupby('chrom').agg({'pos':lambda x: list(x)}).reset_index(drop=True)
        for i in [0,1]:
            try:
                assert np.all([np.all([abs(x-y) <=2000 for y in g.loc[i,'pos']]) for x in g.loc[i,'pos']])
            except: raise AssertionError('Sample: ' + sample + ';Breakpoints are greater than 1kb from each other in the event ',ev)
    elif len(numChrom) == 1: #if there is one chromosome do the same
        sys.exit('Sample: ' + sample + ';There is an intrachromosomal bnd')

#%% 12. create event id and choose columns to split    
#test = svdf_bnd.loc[svdf_bnd.event == 'E774226:contig3',:]
#test2 = svdf_bnd.loc[svdf_bnd.cluster_id == 'E422543',:]
#add event column to svdf (non bnd)
svdf['event'] = [svdf.id[i] + ':' + svdf.cluster_id[i] for i in np.arange(0,svdf.shape[0])] #give event id to non bnd svs
svdf['mates'] = ['NaN']*svdf.shape[0] #give mate column to non bnd svs

#make list of columns that are kept separate between two breakpoints in each event
#columns to drop: mateid,chrom_mate,end
svdf_bnd.columns.values
split_cols = ['id','chrom','pos','quality','filter','precise','imprecise','somatic','germline','unknown_length','consensus_sequence','connection_type',
              'confidence_interval_pos','confidence_interval_end','cluster_id','contig_id','contig_size','reads_used_for_assembly','average_coverage',
              'tumor_bkpt1_depth','tumor_bkpt1_sp_reads','tumor_bkpt1_qual','tumor_bkpt1_high_qual_sp_reads','tumor_bkpt1_high_qual_qual',
              'normal_bkpt1_depth','normal_bkpt1_sp_reads','normal_bkpt1_qual','normal_bkpt1_high_qual_sp_reads','normal_bkpt1_high_qual_qual',
              'tumor_bkpt2_depth','tumor_bkpt2_sp_reads','tumor_bkpt2_qual','tumor_bkpt2_high_qual_sp_reads','tumor_bkpt2_high_qual_qual',
              'normal_bkpt2_depth','normal_bkpt2_sp_reads','normal_bkpt2_qual','normal_bkpt2_high_qual_sp_reads','normal_bkpt2_high_qual_qual',
              'tumor_bkpt1_discordant_reads','normal_bkpt1_discordant_reads','tumor_bkpt2_discordant_reads','normal_bkpt2_discordant_reads']
common_cols_events = ['event','type','length','mates']

#make new column names for values that are kept separate when breakends are combined into an event
split_cols_events_be1 = [[x+sffx for sffx in ['_be1','_be1_alt','_be1_alt2']] for x in split_cols]
split_cols_events_be1 = [y for x in split_cols_events_be1 for y in x]
split_cols_events_be2 = [[x+sffx for sffx in ['_be2','_be2_alt','_be2_alt2']] for x in split_cols]
split_cols_events_be2 = [y for x in split_cols_events_be2 for y in x]

#make list of column names for events df
edf_cols = common_cols_events+split_cols_events_be1+split_cols_events_be2

#%% 13. create event data frame by joining bnd svs into events
#convert to events df
evdf = pd.DataFrame()
rowsList = []

for ix,ev in enumerate(set(svdf_bnd.event)):
    #get breakpoint rows for event
    t = svdf_bnd.loc[svdf_bnd.event==ev,].reset_index(drop=True)
    chrs = set(t.chrom)
    #get ordered list of chromosomes in event (chromosomes is ordered list of all chromosomes)
    order = np.array(chromosomes)[[x in chrs for x in chromosomes]]
    #for each chromosome get breakpoint values
    val = []
    for ch in order:
        c = t.loc[t.chrom==ch].reset_index(drop=True)
        nbp = c.shape[0]
        #if there are two breakpoints for chromosome
        if nbp >= 4:
            #we will drop the sv with the lowest quality score until max 3 are left
            #dropIX = np.arange(0,nbp)[[x == min(c.quality) for x in c.quality]][0] #if there is a tie, the first is dropped
            while nbp >= 4: 
                dropIX = np.arange(0,nbp)[[x == min(c.quality) for x in c.quality]][0] #if there is a tie, the first is dropped
                c = c.loc[[x != dropIX for x in np.arange(0,nbp)]].reset_index(drop=True)
                nbp = c.shape[0]
            v = [[c.loc[i,x] for i in [0,1,2]] for x in split_cols]
            v = [y for x in v for y in x]        
            val = val + v
        elif nbp == 3:
            v = [[c.loc[i,x] for i in [0,1,2]] for x in split_cols]
            v = [y for x in v for y in x]        
            val = val + v
        elif nbp == 2:
            v = [[c.loc[i,x] for i in [0,1]] + ['NaN'] for x in split_cols]
            v = [y for x in v for y in x]        
            val = val + v
        elif nbp == 1:
            v = [[c.loc[0,x],'NaN','NaN'] for x in split_cols]
            v = [y for x in v for y in x]
            val = val + v
        else: sys.exit("event " + ev + " has 0 breakpoints in one breakend "+sample)
    #get common values
    val = [t.loc[0,x] for x in common_cols_events] + val
    #v = [t.loc[0,'event'],t.loc[0,'type'],np.mean(t.qscore),t.loc[0,'filter']]+v
    #add combined breakends as unique event
    r = pd.Series(val,index=edf_cols,name=ix)
    rowsList.append(r)
    
evdf = pd.DataFrame(rowsList)[edf_cols]

#%% 14. combine the events df with non bnd svs
#modify non bnd columns and names and bind with BND events
rename_dict = dict(zip(split_cols,[x+'_be1' for x in split_cols]))
rename_dict['end'] = 'pos_be2' #end coordinate is breakend2 
rename_dict['chrom_mate'] = 'chrom_be2' #mate chrom is breakend2 chrom
svdf = svdf.rename(columns=rename_dict)
#drop mateid
svdf = svdf.drop(['mateid'],axis=1)
#add NaN values for be2 columns
nanFill = ['NaN']*svdf.shape[0]
new_cols_dict = dict(zip([q for p in [[x+y for y in ['_be1_alt','_be1_alt2','_be2','_be2_alt','_be2_alt2']] for x in split_cols] for q in p],[nanFill]*(len(split_cols)*5)))
#don't change pos_be2 or chrom_be2
del new_cols_dict['pos_be2']
del new_cols_dict['chrom_be2']

#add new columns to non bnd structural variants
svdf = pd.concat([svdf,pd.DataFrame(new_cols_dict)],axis=1)

###append bnd snvs to rest
evdf = svdf.append(evdf,ignore_index=True)
###reorder columns
#evdf.columns.values
col_order = ['event','type','mates', #common columns
             'id_be1','id_be1_alt','id_be1_alt2','id_be2','id_be2_alt','id_be2_alt'] + [q for p in [[y+x for y in ['chrom','pos']]
                   for x in ['_be1','_be1_alt','_be1_alt2','_be2','_be2_alt','_be2_alt2']] for q in p] + ['length'] + [q for p in [[y+x for y in ['quality',
              'filter','precise','imprecise','somatic','germline','unknown_length','consensus_sequence','connection_type','confidence_interval_pos',
              'confidence_interval_end','cluster_id','contig_id','contig_size','reads_used_for_assembly','average_coverage','tumor_bkpt1_depth',
              'tumor_bkpt1_sp_reads','tumor_bkpt1_qual','tumor_bkpt1_high_qual_sp_reads','tumor_bkpt1_high_qual_qual','normal_bkpt1_depth',
              'normal_bkpt1_sp_reads','normal_bkpt1_qual','normal_bkpt1_high_qual_sp_reads','normal_bkpt1_high_qual_qual','tumor_bkpt2_depth',
              'tumor_bkpt2_sp_reads','tumor_bkpt2_qual','tumor_bkpt2_high_qual_sp_reads','tumor_bkpt2_high_qual_qual','normal_bkpt2_depth',
              'normal_bkpt2_sp_reads','normal_bkpt2_qual','normal_bkpt2_high_qual_sp_reads','normal_bkpt2_high_qual_qual','tumor_bkpt1_discordant_reads',
              'normal_bkpt1_discordant_reads','tumor_bkpt2_discordant_reads','normal_bkpt2_discordant_reads']] for x in ['_be1','_be1_alt','_be1_alt2',
              '_be2','_be2_alt','_be2_alt2']] for q in p]
evdf = evdf[col_order]
### change data types
for col in ['chrom_be1','chrom_be1_alt','chrom_be1_alt2','chrom_be2','chrom_be2_alt','chrom_be2_alt2']:
    evdf[col] = pd.Categorical(evdf[col], 
                      categories=chromosomes,
                      ordered=True)    
dtype_dict = {'pos_be1':int,'pos_be1_alt':str,'pos_be1_alt2':str,'pos_be2':int,'pos_be2_alt':str,'pos_be2_alt2':str,'length':int}
'''
dtype_dict = {'pos_be1':int,'pos_be1_alt':str,'pos_be1_alt2':str,'pos_be2':int,'pos_be2_alt':str,'pos_be2_alt2':str,
              'precise_be1':bool,'precise_be1_alt':bool,'precise_be1_alt2':bool,'precise_be2':bool,'precise_be2_alt':bool,'precise_be2_alt2':bool,
              'imprecise_be1':bool,'imprecise_be1_alt':bool,'imprecise_be1_alt2':bool,'imprecise_be2':bool,'imprecise_be2_alt':bool,'imprecise_be2_alt2':bool,
              'somatic_be1':bool,'somatic_be1_alt':bool,'somatic_be1_alt2':bool,'somatic_be2':bool,'somatic_be2_alt':bool,'somatic_be2_alt2':bool,
              'germline_be1':bool,'germline_be1_alt':bool,'germline_be1_alt2':bool,'germline_be2':bool,'germline_be2_alt':bool,'germline_be2_alt2':bool,
              'unknown_length_be1':bool,'unknown_length_be1_alt':bool,'unknown_length_be1_alt2':bool,'unknown_length_be2':bool,'unknown_length_be2_alt':bool,'unknown_length_be2_alt2':bool,
              'length':int}
'''
evdf = evdf.astype(dtype_dict)
#order by breakend 1 chrom and pos
evdf = evdf.sort_values(['chrom_be1','pos_be1','chrom_be2','pos_be2']).reset_index(drop=True)

#%% 15. save the events df if desired
#evdf.to_csv(outPath,sep='\t',header=True,index=False)