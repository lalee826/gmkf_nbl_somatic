"""
This script prepares raw VCF files produced by the structural variant callers GridSS, Novobreak, Lumpy, and Manta for downstream analysis
The output is an events table describing SV types with all relevant fields for filtering high confidence calls
"""

import os
import numpy as np
import pandas as pd

workdir = '/rocker-build/gmkf_nbl_somatic/'

chromosomes = ['chr' + str(i) for i in np.arange(1,23)] + ['chrX','chrY']

#%%pull sample metadata
clinData = pd.read_csv(workdir + 'Data/tumor_clinical_data.tsv',sep='\t',header=0)
metadata = pd.read_table(workdir + '/Data/sample_manifest.csv',sep=',',header=0,index_col=None)
#pull KF metadata
pid = [metadata['Kids First Participant ID'].loc[metadata['Kids First Biospecimen ID']==x].reset_index().at[0,'Kids First Participant ID'] for x in clinData['Tumor_Sample_Barcode'].tolist()]
nid = [metadata['Kids First Biospecimen ID'].loc[np.logical_and(metadata['Kids First Participant ID']==x,metadata['sample_type']=='Normal')].reset_index().at[0,'Kids First Biospecimen ID'] for x in pid]
clinData['Normal_Sample_Barcode'] = nid
tumor_normal_conv_table = clinData[['Tumor_Sample_Barcode','Normal_Sample_Barcode']]

#%% get all files in VCF directories
#get files in lumpy_dr
lf = os.listdir(workdir+'Raw Data/Lumpy VCF')
#grab only vcf files
lf = np.array(lf)[np.array([(lambda x: True if x[-3:] == 'vcf' else False)(x) for x in lf])]
#remove any hidden files
lf = np.array(lf)[np.array([(lambda x: False if x.startswith('.') else True)(x) for x in lf])]
#work with raw files only from data in repo
lf = np.array(lf)[np.array([(lambda x: False if 'clean' in x.split('.') else True)(x) for x in lf])]

#get files in manta dir
mf = os.listdir(workdir+'Raw Data/Manta VCF')
mf = np.array(mf)[np.array([(lambda x: True if x[-3:] == 'vcf' else False)(x) for x in mf])]
mf = np.array(mf)[np.array([(lambda x: False if x.startswith('.') else True)(x) for x in mf])]
mf = np.array(mf)[np.array([(lambda x: False if 'clean' in x.split('.') else True)(x) for x in mf])]

#get files in novobreak dirmf = os.listdir(wd+'manta_vcf')
nf = os.listdir(workdir+'Raw Data/Novobreak VCF')
nf = np.array(nf)[np.array([(lambda x: True if x[-3:] == 'vcf' else False)(x) for x in nf])]
nf = np.array(nf)[np.array([(lambda x: False if x.startswith('.') else True)(x) for x in nf])]
nf = np.array(nf)[np.array([(lambda x: False if 'clean' in x.split('.') else True)(x) for x in nf])]

#get files in gridss dir = os.listdir(wd+'gridss_vcf')
gf = os.listdir(workdir+'Raw Data/Gridss VCF')
gf = np.array(gf)[np.array([(lambda x: True if x[-3:] == 'vcf' else False)(x) for x in gf])]
gf = np.array(gf)[np.array([(lambda x: False if x.startswith('.') else True)(x) for x in gf])]
gf = np.array(gf)[np.array([(lambda x: False if 'clean' in x.split('.') else True)(x) for x in gf])]


#%% CLEANING GRIDSS VCFS
'''Function retrieves case_id within VCF file by pulling sample name out of header'''
def getCID(file):
    
    with open(workdir+ 'Raw Data/Gridss VCF/'+file,'r') as infile:
    
        for line in infile.readlines():
            #find the header line
            if line.startswith('#CHROM'): 
                ##the sample barcode can either be the last or second to last element in the header
                ##the sample barcode also can either be the tumor or normal sample barcode
                #check the second to last element first; check last element after
                if line.split('\t')[-2][2] == '_':
                    bc = line.split('\t')[-2].strip('\n')
                    #if tumor, return case_id; if normal, get tumor bc and return case id
                    if bc in tumor_normal_conv_table['Tumor_Sample_Barcode'].tolist():
                        return(clinData.loc[clinData['Tumor_Sample_Barcode']==bc,'case_id'].to_list())
                    elif bc in tumor_normal_conv_table['Normal_Sample_Barcode'].tolist():
                        tbc = tumor_normal_conv_table.loc[tumor_normal_conv_table['Normal_Sample_Barcode']==bc,'Tumor_Sample_Barcode'].tolist()[0]
                        return(clinData.loc[clinData['Tumor_Sample_Barcode']==tbc,'case_id'].to_list())
                elif line.split('\t')[-1][2] == '_':
                    bc = line.split('\t')[-1].strip('\n')
                    #if tumor, return case_id; if normal, get tumor bc and return case id
                    if bc in tumor_normal_conv_table['Tumor_Sample_Barcode'].tolist():
                        return(clinData.loc[clinData['Tumor_Sample_Barcode']==bc,'case_id'].to_list())
                    elif bc in tumor_normal_conv_table['Normal_Sample_Barcode'].tolist():
                        tbc = tumor_normal_conv_table.loc[tumor_normal_conv_table['Normal_Sample_Barcode']==bc,'Tumor_Sample_Barcode'].tolist()[0]
                        return(clinData.loc[clinData['Tumor_Sample_Barcode']==tbc,'case_id'].to_list())
            else: continue
        return([])

#save a new vcf file with headers removed and proper sample id
for file in gf:
    
    #get case_id
    cid = getCID(file)
    if len(cid) == 0:
        print(file)
        continue

    #save file if desired
    '''    
    with open(workdir+'Raw Data/Gridss VCF/'+file,'r') as infile:
        with open(workdir+'Raw Data/Gridss VCF/'+cid[0]+'.gridss.clean.vcf','w') as outfile:        
            for line in infile.readlines():
                
                if line.startswith('##'):
                    continue
                elif line.startswith('#CHROM'):
                    line = line[1:]
                    outfile.write(line)
                elif line.startswith('chr'):
                    linesplit = line.split('\t')
                    if linesplit[0] in chromosomes:
                        line = ('\t').join(linesplit)
                        outfile.write(line)
    '''

#%% convert gridss vcfs to table of calls per sample
#use example file from previous step
df = pd.read_table(workdir+'Raw Data/Gridss VCF/'+'PATBFN'+'.gridss.clean.vcf',sep='\t',header=0)

#%%
#########################################################
#######determine which column is tumor and which is normal
##the sample barcode can either be the last or second to last element in the header
##the sample barcode also can either be the tumor or normal sample barcode
#check the second to last element first; check last element after
if df.columns.values[-2][2] == '_':
    bc = df.columns.values[-2]
    #set col numbers of tumor and normal
    if bc in tumor_normal_conv_table['Tumor_Sample_Barcode'].tolist():
        tumor_col = df.shape[1]-2
        normal_col = df.shape[1]-1
    elif bc in tumor_normal_conv_table['Normal_Sample_Barcode'].tolist():
        tumor_col = df.shape[1]-1
        normal_col = df.shape[1]-2        
elif df.columns.values[-1][2] == '_':
    bc = df.columns.values[-1]
    #set col numbers of tumor and normal
    if bc in tumor_normal_conv_table['Tumor_Sample_Barcode'].tolist():
        tumor_col = df.shape[1]-1
        normal_col = df.shape[1]-2
    elif bc in tumor_normal_conv_table['Normal_Sample_Barcode'].tolist():
        tumor_col = df.shape[1]-2
        normal_col = df.shape[1]-1
else:
    #use metadata file to convert names to tumor sample barcodes
    name = df.columns.values[-1].split('.')[0]
    bc = metadata['Kids First Biospecimen ID'][np.array([x.split('.')[0] == name for x in metadata['name']])].tolist()[0]
    #set col numbers of tumor and normal
    if bc in tumor_normal_conv_table['Tumor_Sample_Barcode'].tolist():
        tumor_col = df.shape[1]-1
        normal_col = df.shape[1]-2
    elif bc in tumor_normal_conv_table['Normal_Sample_Barcode'].tolist():
        tumor_col = df.shape[1]-2
        normal_col = df.shape[1]-1
        
###take columns needed to make variant table
info = df.INFO
formats = df.FORMAT
tumorValues = df.iloc[:,tumor_col]
normalValues = df.iloc[:,normal_col]

###create these to make table of sv calls

info_dicts = [] #metrics for variant call in list of dictionaries
is_imprecise = [] #bool for IMPRECISE metric
tnv = [] #metrics for tumor and normal in list of dictionaries
for i in np.arange(0,len(info)):
    #split info field
    infoSplit = info[i].split(';')
    #split each info field by key,value
    infoSplit = [x.split('=') for x in infoSplit]
    #check if call is imprecise
    isImprecise = np.any([x[0]=='IMPRECISE' for x in infoSplit])
    '''
    #remove imprecise and check if any other fields are len==1
    if np.any([len(y) != 2 for y in list(np.array(infoSplit)[[x[0]!='IMPRECISE' for x in infoSplit]])]):
        print(i)
    '''
    #remove any fields that are not of key=value format
    infoSplit = list(np.array(infoSplit)[[len(x)==2 for x in infoSplit]])
    #create dictionary of fields and values
    #save row to lists
    info_dicts.append(dict(tuple([tuple(x) for x in infoSplit])))
    is_imprecise.append(isImprecise)
    
    #get field keys
    form = formats[i].split(':')
    #get values for tumor and normal
    tv = tumorValues[i].split(':')
    nv = normalValues[i].split(':')
    #create tumor and normal dictionaries
    td = dict(tuple([tuple([form[j],tv[j]]) for j in np.arange(0,len(form))]))
    nd = dict(tuple([tuple([form[j],nv[j]]) for j in np.arange(0,len(form))]))
    tnv.append({'tum':td,'norm':nd})

###extract the values to be included in table
###### PRIMARY METRICS
#pull chroms
chroms = df.CHROM.tolist()
#pull POS
pos = df.POS.tolist()
#pull ALT notation
alt = df.ALT.tolist()
#pull ID
ids = df.ID.tolist() #gridSS call IDs
#pull EVENT ID
events = [info_dicts[i]['EVENT'] for i in np.arange(0,len(info_dicts))]
#pull sv type
#types = [info_dicts[i]['EVENTTYPE'] for i in np.arange(0,len(info_dicts))]
types = [info_dicts[i]['SVTYPE'] for i in np.arange(0,len(info_dicts))]
#pull quality score
qs = df.QUAL.tolist()
#pull filter metric
fil = df.FILTER.tolist()

###### BREAKEND METRICS
beas = [info_dicts[i]['AS'] for i in np.arange(0,len(info_dicts))]
beasrp = [info_dicts[i]['ASRP'] for i in np.arange(0,len(info_dicts))]
beassr = [info_dicts[i]['ASSR'] for i in np.arange(0,len(info_dicts))]
beba = [info_dicts[i]['BA'] for i in np.arange(0,len(info_dicts))]
#bebasrp = [info_dicts[i]['BASRP'] for i in np.arange(0,len(info_dicts))]
#bebassr = [info_dicts[i]['BASSR'] for i in np.arange(0,len(info_dicts))]
bebanrp = [info_dicts[i]['BANRP'] for i in np.arange(0,len(info_dicts))]
bebansr = [info_dicts[i]['BANSR'] for i in np.arange(0,len(info_dicts))]
bebvf = [info_dicts[i]['BVF'] for i in np.arange(0,len(info_dicts))]
beras = [info_dicts[i]['RAS'] for i in np.arange(0,len(info_dicts))]
beref = [info_dicts[i]['REF'] for i in np.arange(0,len(info_dicts))]
berp = [info_dicts[i]['RP'] for i in np.arange(0,len(info_dicts))]
berefpair = [info_dicts[i]['REFPAIR'] for i in np.arange(0,len(info_dicts))]
bevf = [info_dicts[i]['VF'] for i in np.arange(0,len(info_dicts))]
besb = [info_dicts[i]['SB'] for i in np.arange(0,len(info_dicts))]
#besr = [info_dicts[i]['SR'] for i in np.arange(0,len(info_dicts))]
#betaf = [info_dicts[i]['TAF'] for i in np.arange(0,len(info_dicts))]

####### TUMOR METRICS
tasrp = [tnv[i]['tum']['ASRP'] for i in np.arange(0,len(info_dicts))]
tassr = [tnv[i]['tum']['ASSR'] for i in np.arange(0,len(info_dicts))]
tbanrp = [tnv[i]['tum']['BANRP'] for i in np.arange(0,len(info_dicts))]
tbansr = [tnv[i]['tum']['BANSR'] for i in np.arange(0,len(info_dicts))]
#tbasrp = [tnv[i]['tum']['BASRP'] for i in np.arange(0,len(info_dicts))]
#tbassr = [tnv[i]['tum']['BASSR'] for i in np.arange(0,len(info_dicts))]
tbvf = [tnv[i]['tum']['BVF'] for i in np.arange(0,len(info_dicts))]
tref = [tnv[i]['tum']['REF'] for i in np.arange(0,len(info_dicts))]
trp = [tnv[i]['tum']['RP'] for i in np.arange(0,len(info_dicts))]
trefpair = [tnv[i]['tum']['REFPAIR'] for i in np.arange(0,len(info_dicts))]
tvf = [tnv[i]['tum']['VF'] for i in np.arange(0,len(info_dicts))]

####### NORMAL METRICS
nasrp = [tnv[i]['norm']['ASRP'] for i in np.arange(0,len(info_dicts))]
nassr = [tnv[i]['norm']['ASSR'] for i in np.arange(0,len(info_dicts))]
nbanrp = [tnv[i]['norm']['BANRP'] for i in np.arange(0,len(info_dicts))]
nbansr = [tnv[i]['norm']['BANSR'] for i in np.arange(0,len(info_dicts))]
#nbasrp = [tnv[i]['norm']['BASRP'] for i in np.arange(0,len(info_dicts))]
#nbassr = [tnv[i]['norm']['BASSR'] for i in np.arange(0,len(info_dicts))]
nbvf = [tnv[i]['norm']['BVF'] for i in np.arange(0,len(info_dicts))]
nref = [tnv[i]['norm']['REF'] for i in np.arange(0,len(info_dicts))]
nrp = [tnv[i]['norm']['RP'] for i in np.arange(0,len(info_dicts))]
nrefpair = [tnv[i]['norm']['REFPAIR'] for i in np.arange(0,len(info_dicts))]
nvf = [tnv[i]['norm']['VF'] for i in np.arange(0,len(info_dicts))]

#create a data frame to output
svdf = pd.DataFrame({'event':events,'id':ids,'type':types,'qscore':qs,'filter':fil,'chrom':chroms,'pos':pos,'alt':alt,'imprecise':is_imprecise, #main 
                     'be_support_assembly_tot':beas,'be_support_assembly_incorp_rp':beasrp,'be_support_assembly_sr':beassr,'be_support_assembly_local':beba,'be_support_assembly_re':beras, #assembly support on breakend
                     'be_support_fragments_bp':bebvf,'be_support_fragments_alt':bevf,                   #fragments supporting
                     'be_mapacross_tot':beref,'be_mapacross_alt':berp,'be_mapacross_ref':berefpair,     #map across junction metrics
                     'be_nonsupport_contig_rp':bebanrp,'be_nonsupport_contig_sr':bebansr,               #metrics for reads not supporting
                     'be_strandbias':besb,#'be_taf':betaf,                                               #misc
                     'tum_support_assembly_incorp_rp':tasrp,'tum_support_assembly_sr':tassr,            #breakend support in tumor
                     'tum_support_fragments_bp':tbvf,'tum_support_fragments_alt':tvf,                   #fragments supporting tumor
                     'tum_mapacross_tot':tref,'tum_mapacross_alt':trp,'tum_mapacross_ref':trefpair,     #map across junction metrics tumor
                     'tum_nonsupport_contig_rp':tbanrp,'tum_nonsupport_contig_sr':tbansr,               #metrics for reads not supporting in tumor
                     'norm_support_assembly_incorp_rp':nasrp,'norm_support_assembly_sr':nassr,          #breakend support in NORMAL
                     'norm_support_fragments_bp':nbvf,'norm_support_fragments_alt':nvf,                 #fragments supporting NORMAL
                     'norm_mapacross_tot':nref,'norm_mapacross_alt':nrp,'norm_mapacross_ref':nrefpair,  #map across junction metrics NORMAL
                     'norm_nonsupport_contig_rp':nbanrp,'norm_nonsupport_contig_sr':nbansr             #metrics for reads not supporting in NORMAL
                     })
svdf.chrom = pd.Categorical(svdf.chrom, 
                      categories=['chr'+str(x) for x in list(range(1,23))+['X','Y']],
                      ordered=True)
#svdf = svdf.sort_values(by=['event','chrom','pos'])
#make list of columns that are kept separte between two breakpoints in each event
split_cols = ['id','chrom','pos','alt','imprecise','be_support_assembly_tot', 'be_support_assembly_incorp_rp',
   'be_support_assembly_sr', 'be_support_assembly_local',
   'be_support_assembly_re', 'be_support_fragments_bp',
   'be_support_fragments_alt', 'be_mapacross_tot', 'be_mapacross_alt',
   'be_mapacross_ref', 'be_nonsupport_contig_rp',
   'be_nonsupport_contig_sr', 'be_strandbias',# 'be_taf',
   'tum_support_assembly_incorp_rp', 'tum_support_assembly_sr',
   'tum_support_fragments_bp', 'tum_support_fragments_alt',
   'tum_mapacross_tot', 'tum_mapacross_alt', 'tum_mapacross_ref',
   'tum_nonsupport_contig_rp', 'tum_nonsupport_contig_sr',
   'norm_support_assembly_incorp_rp', 'norm_support_assembly_sr',
   'norm_support_fragments_bp', 'norm_support_fragments_alt',
   'norm_mapacross_tot', 'norm_mapacross_alt', 'norm_mapacross_ref',
   'norm_nonsupport_contig_rp', 'norm_nonsupport_contig_sr']
#make new column names for values that are kept separate when breakends are combined into an event
split_cols_events = [[x+sffx for sffx in ['_be1','_be2']] for x in split_cols]
split_cols_events = [y for x in split_cols_events for y in x]
#make list of column names for events df
edf_cols = ['event','type','qscore','filter']+split_cols_events


############# keep only variants that PASS filter
svdf_pass = svdf.loc[svdf['filter'] == 'PASS'].reset_index(drop=True)
##keep track of non-mated events
nonMateEvents = []

#%%create dictionary of event names to row indices for quick access
evList = svdf_pass.event
evIX = svdf_pass.index.values
evIXDict = {}
addEvDict = np.vectorize(lambda x: evIXDict.update({evList[x]:[evIX[x]]}) if evList[x] not in evIXDict else evIXDict[evList[x]].append(evIX[x]))
addEvDict(np.arange(len(evIX)))
#save each row in list
rowsList = []
i = 0
for ix,ev in enumerate(set(svdf_pass.event)):
    '''
    ev = 'gridss0_22103'
    '''
    #get breakpoint rows for event
    evIXS = evIXDict[ev]
    t = svdf_pass.loc[evIXS,].reset_index(drop=True)
    #skip if bnd is unmated
    if t.shape[0] == 1:
        #print(ev,",",t['type'][0])
        nonMateEvents.append(ev)
        continue
    #sort by chrom and pos
    t = t.sort_values(['chrom','pos'])
    #get values for break ends 1 and 2 for split columns
    v = [[t.loc[i,x] for i in [0,1]] for x in split_cols]
    v = [y for x in v for y in x]
    #get common event,type,combined qscore,and common filter value
    v = [t.loc[0,'event'],t.loc[0,'type'],np.mean(t.qscore),t.loc[0,'filter']]+v
    #add combined breakends as unique event
    r = pd.Series(v,index=edf_cols,name=i)
    i += 1
    #save row to list
    rowsList.append(r)
    #evdf = evdf.append(r)[edf_cols]
#%%#convert to events df
evdf = pd.DataFrame(rowsList)[edf_cols]
evdf.pos_be1 = evdf.pos_be1.astype(int)
evdf.pos_be2 = evdf.pos_be2.astype(int)
n = evdf.shape[0]
#separate BND from rest
bndIX = np.arange(n)[[x!=y for x,y in zip(evdf.chrom_be1,evdf.chrom_be2)]]
nonIX = np.arange(n)[[x==y for x,y in zip(evdf.chrom_be1,evdf.chrom_be2)]]
nonIXDF = evdf.loc[nonIX][['alt_be1','alt_be2','pos_be1','pos_be2']]
#call duplications
dupIX = nonIX[[len(x.split(']'))==3 and len(y.split('['))==3 for x,y in zip(nonIXDF.alt_be1,nonIXDF.alt_be2)]]
#call inversions
invIX = nonIX[[np.logical_or(len(x.split('['))==3 and len(y.split('['))==3,len(x.split(']'))==3 and len(y.split(']'))==3) for x,y in zip(nonIXDF.alt_be1,nonIXDF.alt_be2)]]
#call deletions/insertions
delInsIX = nonIX[[len(x.split('['))==3 and len(y.split(']'))==3 for x,y in zip(nonIXDF.alt_be1,nonIXDF.alt_be2)]]
delInsDF = nonIXDF.loc[delInsIX]
#insertions have difference of 1 between positions
delIx = delInsIX[[y-x > 32 for x,y in zip(delInsDF.pos_be1,delInsDF.pos_be2)]]
insIx = delInsIX[[y-x <= 32 for x,y in zip(delInsDF.pos_be1,delInsDF.pos_be2)]]
#set types by index
evdf.loc[dupIX,'type'] = 'DUP'
evdf.loc[invIX,'type'] = 'INV'
evdf.loc[delIx,'type'] = 'DEL'
evdf.loc[insIx,'type'] = 'INS'
sizes = []
#add size column if not BND
for ix,row in evdf.iterrows():
    t = row['type']
    if t in ['BND','INV']: sizes.append('NA')
    else:
        sizes.append(row['pos_be2']-row['pos_be1'])
evdf['size'] = sizes

evdfDel = evdf.loc[evdf.type=='DEL']
evdfIns = evdf.loc[evdf.type=='INS']

#%% save gridss output
evdf.chrom_be1 = pd.Categorical(evdf.chrom_be1, 
                      categories=['chr'+str(x) for x in list(range(1,23))+['X','Y']],
                      ordered=True)
evdf.chrom_be2 = pd.Categorical(evdf.chrom_be2, 
                      categories=['chr'+str(x) for x in list(range(1,23))+['X','Y']],
                      ordered=True)
evdf = evdf.sort_values(['chrom_be1','pos_be1']).reset_index(drop=True)

##save output if desired
#evdf.to_csv('',sep='\t',header=True,index=False)

#%% CLEANING LUMPY VCFS 
#remove header lines for each sample and collapse to one data table
def getCID(inFile):
    
    with open(workdir+'Raw Data/Lumpy VCF/'+inFile,'r') as infile:
    
        for line in infile.readlines():
                
            if line.startswith('#CHROM'):
                tbc = line.split('\t')[-2]
                cid = clinData.loc[clinData['Tumor_Sample_Barcode']==tbc,'case_id'].to_list()
                if len(cid) > 0:
                    cid = cid[0]
                return(cid)
            else: continue


for file in lf:

    #get case_id
    cid = getCID(file)
    if len(cid) == 0:
        print(file)
        continue
    
    #save clean vcf if desired
    '''
    with open(workdir+'Raw Data/Lumpy VCF/'+file,'r') as infile:
        with open(workdir+'Raw Data/Lumpy VCF/'+cid+'.lumpy.clean.vcf','w') as outfile:        
            for line in infile.readlines():
                
                if line.startswith('##'):
                    continue
                elif line.startswith('#CHROM'):
                    line = line[1:]
                    outfile.write(line)
                elif line.startswith('chr'):
                    linesplit = line.split('\t')
                    if linesplit[0] in chromosomes:
                        line = ('\t').join(linesplit)
                        outfile.write(line)
    '''

#%% convert clean lumpy vcf to table of calls per sample
#use example file in repo
df = pd.read_table((workdir+'Raw Data/Lumpy VCF/PATBFN.lumpy.clean.vcf'),sep='\t',header=0)

info = df['INFO'].tolist()
info = [x.split(';') for x in info]

svtype = [x[0].split('=')[1] for x in info] #get sv types

#get sv lengths
lens = [np.arange(len(x))[[y.startswith('SVLEN') for y in x]] for x in info]
lens = [(lambda x: int(x) if len(x) != 0 else None)(x) for x in lens]
lens = [(info[i][lens[i]]).split('=')[1] if lens[i] != None else None for i in np.arange(len(lens))]
lens = [abs(int(x)) if x != None else None for x in lens]

#get sv ends
ends = [np.arange(len(x))[[y.startswith('END') for y in x]] for x in info]
ends = [(lambda x: int(x) if len(x) != 0 else None)(x) for x in ends]
ends = [(info[i][ends[i]]).split('=')[1] if ends[i] != None else None for i in np.arange(len(ends))]
ends = [int(x) if x != None else None for x in ends]

#get event
event = [np.arange(len(x))[[y.startswith('EVENT') for y in x]] for x in info]
event = [(lambda x: int(x) if len(x) != 0 else None)(x) for x in event]
event = [(info[i][event[i]]).split('=')[1] if event[i] != None else None for i in np.arange(len(event))]
event = [int(x) if x != None else None for x in event]

#get mate_id
mid = [np.arange(len(x))[[y.startswith('MATEID') for y in x]] for x in info]
mid = [(lambda x: int(x) if len(x) != 0 else None)(x) for x in mid]
mid = [(info[i][mid[i]]).split('=')[1] if mid[i] != None else None for i in np.arange(len(mid))]
mid = [str(x) if x != None else None for x in mid]

#get sv imprecise bool
imp = ['IMPRECISE' in x for x in info]

#get secondary bool
sec = ['SECONDARY' in x for x in info]

#get supporting count for tumor
tum = df.iloc[:,9].tolist()
tum = [x.split(':') for x in tum]
tumor_total_support = [int(x[1]) for x in tum]
tumor_paired_end = [int(x[2]) for x in tum]
tumor_split = [int(x[3]) for x in tum]

#get supporting count for normal
norm = df.iloc[:,10].tolist()
norm = [x.split(':') for x in norm]
normal_total_support = [int(x[1]) for x in norm]
normal_paired_end = [int(x[2]) for x in norm]
normal_split = [int(x[3]) for x in norm]

tdf = pd.DataFrame({'Chromosome':df['CHROM'],
              'Start':df['POS'],
              'End':ends,
              'Length':lens,
              'Type':svtype,
              'Event':event,
              'ID':df['ID'],
              'Alt':df['ALT'],
              'Mate_ID':mid,
              'Is_imprecise':imp,
              'Is_secondary':sec,
              'Tumor_total_support':tumor_total_support,
              'Tumor_paired_end':tumor_paired_end,
              'Tumor_split_reads':tumor_split,
              'Normal_total_support':normal_total_support,
              'Normal_paired_end':normal_paired_end,
              'Normal_split_reads':normal_split
              })

#%% filter for high confidence somatic calls
#filter calls that have sufficient support in tumor and at least one paired end read
#total support >= 10 or paired end support >= 5
keep_tum = np.logical_and(np.logical_or(tdf['Tumor_total_support'] >= 10,tdf['Tumor_paired_end'] >= 5),tdf['Tumor_paired_end'] > 0)

#filter calls that have no support in normal
keep_norm = np.logical_or(tdf['Normal_total_support'] == 0,
                          np.logical_and(tdf['Normal_total_support'] <= 10,
                                         np.logical_and(tdf['Normal_total_support']/tdf['Tumor_total_support'] < 0.1,
                                         tdf['Tumor_total_support'] > 20)))

#combine and filter
tdff = tdf.ix[np.logical_and(keep_tum,keep_norm)]

#save the data frame if desired
#tdff.to_csv((),sep='\t',index=False)


#%% CLEAN MANTA VCFS
#%% remove header lines for each sample and collapse to one data table
def getCID(inFile):
    
    with open(workdir+'Raw Data/Manta VCF/'+inFile,'r') as infile:
    
        for line in infile.readlines():
                
            if line.startswith('#CHROM'):
                tbc = line.split('\t')[-2]
                cid = clinData.loc[clinData['Tumor_Sample_Barcode']==tbc,'case_id'].to_list()
                if len(cid) > 0:
                    cid = cid[0]
                    return(cid)
                else:
                    tbc = line.split('\t')[-1].strip('\n')
                    cid = clinData.loc[clinData['Tumor_Sample_Barcode']==tbc,'case_id'].to_list()
                    if len(cid) > 0:
                        cid = cid[0]
                    return(cid)
            else: continue

for file in mf:

    #get case_id
    cid = getCID(file)
    if len(cid) == 0:
        print(file)
        continue
    
    #save cleaned vcf if desired
    '''
    with open(workdir+'/Raw Data/Manta VCF/'+file,'r') as infile:
        with open(workdir+'/Raw Data/Manta VCF/'+cid+'.manta.clean.vcf','w') as outfile:        
            for line in infile.readlines():
                
                if line.startswith('##'):
                    continue
                elif line.startswith('#CHROM'):
                    line = line[1:]
                    outfile.write(line)
                elif line.startswith('chr'):
                    linesplit = line.split('\t')
                    if linesplit[0] in chromosomes:
                        line = ('\t').join(linesplit)
                        outfile.write(line)
    '''


df = pd.read_table(workdir+'Raw Data/Manta VCF/PATBFN.manta.clean.vcf',sep='\t',header=0)
    
#get chromosomes
chrom = df['CHROM']

#get start positions
starts = df['POS']

#get var type
vtypes = [x.strip('<') for x in df['ALT']]
vtypes = [x.strip('>') for x in vtypes]
vtypes = ['DEL' if len(x) < 3 else x for x in vtypes]

#get var id
vid = df['ID']

#split info field
info = [np.array(x.split(';')) for x in df['INFO']]
#remove csq
infokeep = [np.arange(len(x))[np.invert(np.array([y.startswith('CSQ') for y in x]))] for x in info]
info = [(info[i][infokeep[i]]).tolist() for i in np.arange(len(info))]

#get end pos
ends = [np.arange(len(x))[[y.startswith('END') for y in x]] for x in info]
ends = [(lambda x: int(x) if len(x) != 0 else None)(x) for x in ends]
ends = [(info[i][ends[i]]).split('=')[1] if ends[i] != None else None for i in np.arange(len(ends))]
ends = [abs(int(x)) if x != None else None for x in ends]

#get length
svlen = [np.arange(len(x))[[y.startswith('SVLEN') for y in x]] for x in info]
svlen = [(lambda x: int(x) if len(x) != 0 else None)(x) for x in svlen]
svlen = [(info[i][svlen[i]]).split('=')[1] if svlen[i] != None else None for i in np.arange(len(svlen))]
svlen = [abs(int(x)) if x != None else None for x in svlen]


#get mateid
mid = [np.arange(len(x))[[y.startswith('MATEID') for y in x]] for x in info]
mid = [(lambda x: int(x) if len(x) != 0 else None)(x) for x in mid]
mid = [(info[i][mid[i]]).split('=')[1] if mid[i] != None else None for i in np.arange(len(mid))]

#get sv imprecise bool
imp = ['IMPRECISE' in x for x in info]

#get secondary bool
som = ['SOMATIC' in x for x in info]

tumor_paired_end = []
tumor_split = []
normal_paired_end = []
normal_split = []

#get supporting paired and split reads in tumor
if df.columns.values[-2] in set(clinData['Tumor_Sample_Barcode']):
    tum_vals = df.iloc[:,-2]
    norm_vals = df.iloc[:,-1]
elif df.columns.values[-1] in set(clinData['Tumor_Sample_Barcode']):
    tum_vals = df.iloc[:,-1]
    norm_vals = df.iloc[:,-2]
for i in np.arange(df.shape[0]):

    form = df.loc[i,'FORMAT'].split(':')

    if 'PR' not in form:
        tumor_paired_end.append(0)
        normal_paired_end.append(0)

    if 'SR' not in form:
        tumor_split.append(0)
        normal_split.append(0)

    for j in np.arange(len(form)):

        if form[j] == 'PR':
            tumor_paired_end.append(int(tum_vals[i].split(':')[j].split(',')[1]))
            normal_paired_end.append(int(norm_vals[i].split(':')[j].split(',')[1]))
        elif form[j] == 'SR':
            tumor_split.append(int(tum_vals[i].split(':')[j].split(',')[1]))
            normal_split.append(int(norm_vals[i].split(':')[j].split(',')[1]))

dft = pd.DataFrame({'Chromosome':chrom,
                    'Start':starts,
                    'End':ends,
                    'Length':svlen,
                    'Type':vtypes,
                    'ID':vid,
                    'Mate_ID':mid,
                    'Is_imprecise':imp,
                    'Is_somatic':som,
                    'Tumor_paired_end':tumor_paired_end,
                    'Tumor_split_reads':tumor_split,
                    'Normal_paired_end':normal_paired_end,
                    'Normal_split_reads':normal_split
                        })

#save clean manta VCF if desired
#dft.to_csv(''),sep='\t',index=False)

#%% Clean Novobreak VCFs
#we use a manifest file to export sample data for novobreak
nf_man = pd.read_csv(workdir+'Raw Data/Novobreak VCF/novobreak_manifest.csv')

#novobreak has additional data without headers that should be added in
otherHeaders = ('\t').join(('cluster_id, contig_id, contig_size, reads_used_for_assembly, average_coverage, tumor_bkpt1_depth, tumor_bkpt1_sp_reads,'
                     'tumor_bkpt1_qual, tumor_bkpt1_high_qual_sp_reads, tumor_bkpt1_high_qual_qual, normal_bkpt1_depth, normal_bkpt1_sp_reads,'
                     'normal_bkpt1_qual, normal_bkpt1_high_qual_sp_reads, normal_bkpt1_high_qual_qual, tumor_bkpt2_depth, tumor_bkpt2_sp_reads,'
                     'tumor_bkpt2_qual, tumor_bkpt2_high_qual_sp_reads, tumor_bkpt2_high_qual_qual, normal_bkpt2_depth, normal_bkpt2_sp_reads,'
                     'normal_bkpt2_qual, normal_bkpt2_high_qual_sp_reads, normal_bkpt2_high_qual_qual, tumor_bkpt1_discordant_reads,'
                     'normal_bkpt1_discordant_reads, tumor_bkpt2_discordant_reads, normal_bkpt2_discordant_reads').strip('\n').replace(' ','').split(','))

for file in nf:

    #get case_id
    cid = nf_man.loc[nf_man['name'] == file,'case_id'].values[0][:6]
    
    #save cleaned vcf if desired
    '''
    with open(workdir+'Raw Data/Novobreak VCF/'+file,'r') as infile:
        with open(workdir+'Raw Data/Novobreak VCF/'+cid+'.novobreak.clean.vcf','w') as outfile:        
            for line in infile.readlines():
                
                if line.startswith('##'):
                    continue
                elif line.startswith('#CHROM'):
                    line = line[1:]
                    #add rest of headers
                    line = (line + '\t' + otherHeaders).replace('\n','') + '\n'
                    outfile.write(line)
                elif line.startswith('chr'):
                    linesplit = line.split('\t')
                    if linesplit[0] in chromosomes:
                        #linesplit = linesplit[:10]
                        line = ('\t').join(linesplit)
                        outfile.write(line)
    '''

#%% remake and clean data
#create table of calls for example in repo 
df = pd.read_table(workdir+'Raw Data/Novobreak VCF/PATBFN.novobreak.clean.vcf',sep='\t',header=0)

#get chromosomes
chrom = df['CHROM']

#get start positions
starts = df['POS']

#get var type
vtypes = [x.strip('<') for x in df['ALT']]
vtypes = [x.strip('>') for x in vtypes]
vtypes = ['BND' if x == 'TRA' else x for x in vtypes]

#split info field
info = [np.array(x.split(';')) for x in df['INFO']]
#remove consensus
infokeep = [np.arange(len(x))[np.invert(np.array([y.startswith('CONSENSUS') for y in x]))] for x in info]
info = [(info[i][infokeep[i]]).tolist() for i in np.arange(len(info))]

#get end pos
ends = [np.arange(len(x))[[y.startswith('END') for y in x]] for x in info]
ends = [(lambda x: int(x) if len(x) != 0 else None)(x) for x in ends]
ends = [(info[i][ends[i]]).split('=')[1] if ends[i] != None else None for i in np.arange(len(ends))]
ends = [abs(int(x)) if x != None else None for x in ends]

#get length
svlen = [np.arange(len(x))[[y.startswith('SVLEN') for y in x]] for x in info]
svlen = [(lambda x: int(x) if len(x) != 0 else None)(x) for x in svlen]
svlen = [(info[i][svlen[i]]).split('=')[1] if svlen[i] != None else None for i in np.arange(len(svlen))]
svlen = [abs(int(x)) if x != None else None for x in svlen]

#get mate chromosome
mid = [np.arange(len(x))[[y.startswith('CHR2') for y in x]] for x in info]
mid = [(lambda x: int(x) if len(x) != 0 else None)(x) for x in mid]
mid = [(info[i][mid[i]]).split('=')[1] if mid[i] != None else None for i in np.arange(len(mid))]

#get sv imprecise bool
imp = ['IMPRECISE' in x for x in info]

#get somatic bool
som = ['SOMATIC' in x for x in info]

#give ID
sample = file.split('_')[0][:6]
vid = [sample+':novobreak:'+vtypes[i]+':'+str(i) for i in np.arange(df.shape[0])]

dft = pd.DataFrame({'Chromosome':chrom,
                    'Start':starts,
                    'End':ends,
                    'Length':svlen,
                    'Type':vtypes,
                    'ID':vid,
                    'Mate_chr':mid,
                    'Is_imprecise':imp,
                    'Is_somatic':som,
                        })        

#save output if desired
#dft.to_csv((),sep='\t',index=False)