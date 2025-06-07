"""
Description: Parse SV calls from LUMPY to create a table of SV events. This script will:
    
    - Join Breakpoints into events with multiple ends
    - Tabulate significant metrics
    - Provide a tabular view of final somatic filtered SV calls made by Lumpy
    - The final events tables can be filtered based on the metrics/counts created in table features
"""

#%% import modules and data
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
df = pd.read_table('/rocker-build/gmkf_nbl_somatic/Raw Data/Lumpy VCF/PATBFN.lumpy.clean.vcf',sep='\t',header=0)

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

#check columns in vcf
df.columns.values

###extract data from info and tumor/normal columns
info = df.INFO
formats = df.FORMAT
tumorValues = df.iloc[:,tumor_col]
normalValues = df.iloc[:,normal_col]

###create these to make table of sv calls
info_dicts = [] #metrics for variant call in list of dictionaries
is_imprecise = [] #bool for IMPRECISE metric
is_secondary = [] #bool for SECONDARY metric
tnv = [] #metrics for tumor and normal in list of dictionaries

for i in np.arange(0,len(info)):
    #split info field
    infoSplit = info[i].split(';')
    #split each info field by key,value
    infoSplit = [x.split('=') for x in infoSplit]
    #check if call is imprecise
    isImprecise = np.any([x[0]=='IMPRECISE' for x in infoSplit])
    #check if call is secondary
    isSecondary = np.any([x[0]=='SECONDARY' for x in infoSplit])
    
    '''
    #remove imprecise and check if any other fields are len==1
    if np.any([len(y) != 2 for y in list(np.array(infoSplit)[[x[0] not in ['IMPRECISE','SECONDARY'] for x in infoSplit]])]):
        print(i)
    '''
    
    #remove any fields that are not of key=value format
    infoSplit = list(np.array(infoSplit)[[len(x)==2 for x in infoSplit]])
    #create dictionary of fields and values
    #save row to lists
    info_dicts.append(dict(tuple([tuple(x) for x in infoSplit])))
    is_imprecise.append(isImprecise)
    is_secondary.append(isSecondary)
    
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
#pull end pos
ends = []
for i in np.arange(0,len(info_dicts)):
    try:
        ends.append(info_dicts[i]['END'])
    except KeyError:
        ends.append('NaN')
#pull ID
ids = df.ID.tolist() #lumpy call IDs
#pull EVENT ID
events = []
for i in np.arange(0,len(info_dicts)):
    try:
        events.append(info_dicts[i]['EVENT'])
    except KeyError:
        events.append('NaN')
#pull sv type
types = [info_dicts[i]['SVTYPE'] for i in np.arange(0,len(info_dicts))]
#pull mate ID
mate = []
for i in np.arange(0,len(info_dicts)):
    try:
        mate.append(info_dicts[i]['MATEID'])
    except KeyError:
        mate.append('NaN')
#pull sv length
svlen = []
for i in np.arange(0,len(info_dicts)):
    try:
        svlen.append(abs(int(info_dicts[i]['SVLEN'])))
    except KeyError:
        svlen.append('NaN')
#pull quality score
qs = df.QUAL.tolist()
#pull filter metric
fil = df.FILTER.tolist()

#get all possible info tags
t = [list(info_dicts[i].keys()) for i in np.arange(0,len(info_dicts))]
info_tags = set([item for sublist in t for item in sublist])

###### BREAKEND METRICS
#PE: Number of paired end reads supporting variant across all samples
#SR: Number of split reads supporting variant across all samples
#SU: Number of pieces of evidence supporting the variant across all samples
#STRANDS: Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)
#CIPOS: Confidence interval around POS for imprecise variants
#CIPOS95: Confidence interval (95%) around POS for imprecise variants
#CIEND: Confidence interval around END for imprecise variants
#CIEND95: Confidence interval (95%) around END for imprecise variantsdead
pe = [info_dicts[i]['PE'] for i in np.arange(0,len(info_dicts))]
sr = [info_dicts[i]['SR'] for i in np.arange(0,len(info_dicts))]
su = [info_dicts[i]['SU'] for i in np.arange(0,len(info_dicts))]
cipos = [info_dicts[i]['CIPOS'] for i in np.arange(0,len(info_dicts))]
ciend = [info_dicts[i]['CIEND'] for i in np.arange(0,len(info_dicts))]
cipos95 = [info_dicts[i]['CIPOS95'] for i in np.arange(0,len(info_dicts))]
ciend95 = [info_dicts[i]['CIEND95'] for i in np.arange(0,len(info_dicts))]

#get all possible tags for tumor/normal
tgt = [tnv[i]['tum']['GT'] for i in np.arange(0,len(info_dicts))]
tsu = [tnv[i]['tum']['SU'] for i in np.arange(0,len(info_dicts))]
tpe = [tnv[i]['tum']['PE'] for i in np.arange(0,len(info_dicts))]
tsr = [tnv[i]['tum']['SR'] for i in np.arange(0,len(info_dicts))]

ngt = [tnv[i]['norm']['GT'] for i in np.arange(0,len(info_dicts))]
nsu = [tnv[i]['norm']['SU'] for i in np.arange(0,len(info_dicts))]
npe = [tnv[i]['norm']['PE'] for i in np.arange(0,len(info_dicts))]
nsr = [tnv[i]['norm']['SR'] for i in np.arange(0,len(info_dicts))]

#create a data frame to output
svdf = pd.DataFrame({'event':events,'id':ids,'mateid':mate,'type':types,'chrom':chroms,'pos':pos,'end':ends,'length':svlen, #primary metrics
                     'imprecise':is_imprecise,'secondary':is_secondary, #primary metrics
                     'junction_support_pe':pe,'junction_support_sr':sr,'junction_support_su':su, #junction metrics
                     'confidence_interval_pos':cipos,'confidence_interval_pos_95':cipos95,'confidence_interval_end':ciend,'confidence_interval_end_95':ciend95,
                     'tum_support_pe':tpe,'tum_support_sr':tsr,'tum_support_su':tsu,
                     'norm_support_pe':npe,'norm_support_sr':nsr,'norm_support_su':nsu})
#sort by chromosome
svdf.chrom = pd.Categorical(svdf.chrom, 
                      categories=['chr'+str(x) for x in list(range(1,23))+['X','Y']],
                      ordered=True)

#### For BND events join mate pairs, discard those without mates
#split off BND from rest
svdf_bnd = svdf.loc[svdf.type=='BND',:].reset_index(drop=True)
svdf = svdf.loc[svdf.type!='BND',:].reset_index(drop=True)
#discard rows whose mate pair is not in data frame
allID = svdf_bnd.id.tolist()
svdf_bnd = svdf_bnd.loc[[i in allID for i in svdf_bnd.mateid],:].reset_index(drop=True)

#make list of columns that are kept separate between two breakpoints in each event
split_cols = ['id','chrom','pos','imprecise','secondary',
              'confidence_interval_pos','confidence_interval_pos_95','confidence_interval_end','confidence_interval_end_95']
#make new column names for values that are kept separate when breakends are combined into an event
split_cols_events = [[x+sffx for sffx in ['_be1','_be2']] for x in split_cols]
split_cols_events = [y for x in split_cols_events for y in x]
#make list of column names for events df
common_cols_events = ['event','type','length','junction_support_pe','junction_support_sr','junction_support_su',
            'tum_support_pe','tum_support_sr','tum_support_su',
            'norm_support_pe','norm_support_sr','norm_support_su']
edf_cols = common_cols_events+split_cols_events


#convert to events df
evdf = pd.DataFrame()
for ix,ev in enumerate(set(svdf_bnd.event)):
    #get breakpoint rows for event
    t = svdf_bnd.loc[svdf_bnd.event==ev,].reset_index(drop=True)
    #sort by chrom and pos
    t = t.sort_values(['chrom','pos'])
    #get values for break ends 1 and 2 for split columns
    v = [[t.loc[i,x] for i in [0,1]] for x in split_cols]
    v = [y for x in v for y in x]
    #get common values
    v = [t.loc[0,x] for x in common_cols_events] + v
    #v = [t.loc[0,'event'],t.loc[0,'type'],np.mean(t.qscore),t.loc[0,'filter']]+v
    #add combined breakends as unique event
    r = pd.Series(v,index=edf_cols,name=ix)
    evdf = evdf.append(r)[edf_cols]

evdf.pos_be1 = evdf.pos_be1.astype(int)
evdf.pos_be2 = evdf.pos_be2.astype(int)
evdf.secondary_be1 = evdf.secondary_be1.astype(bool)
evdf.secondary_be2 = evdf.secondary_be2.astype(bool)
evdf.imprecise_be1 = evdf.imprecise_be1.astype(bool)
evdf.imprecise_be2 = evdf.imprecise_be2.astype(bool)

#modify non bnd columns and names and bind with BND events
svdf = svdf.rename(columns={'id':'id_be1','chrom':'chrom_be1','pos':'pos_be1','imprecise':'imprecise_be1','secondary':'secondary_be1',
                     'confidence_interval_pos':'confidence_interval_pos_be1','confidence_interval_pos_95':'confidence_interval_pos_95_be1',
                     'confidence_interval_end':'confidence_interval_end_be1','confidence_interval_end_95':'confidence_interval_end_95_be1',
                     'mateid':'id_be2','end':'pos_be2'})
#chromosome on be2 is same as be1
svdf['chrom_be2'] = svdf.chrom_be1
#add NaN values for be2 columns
nanFill = ['NaN']*svdf.shape[0]
new_cols = {'imprecise_be2':nanFill,'secondary_be2':nanFill,
              'confidence_interval_pos_be2':nanFill,'confidence_interval_pos_95_be2':nanFill,
              'confidence_interval_end_be2':nanFill,'confidence_interval_end_95_be2':nanFill}
svdf = pd.concat([svdf,pd.DataFrame(new_cols)],axis=1)
#give the non BND svs an event id (just take the id)
svdf['event'] = svdf['id_be1']
#order by breakend 1 chrom and pos
evdf = svdf.append(evdf,ignore_index=True)
evdf = evdf.sort_values(['chrom_be1','pos_be1']).reset_index(drop=True)
#reorder columns
evdf.columns.values
evdf = evdf[['event', 'id_be1', 'id_be2', 'type', 'chrom_be1', 'pos_be1', 'chrom_be2', 'pos_be2', 'length', 'imprecise_be1', 'imprecise_be2',
     'secondary_be1', 'secondary_be2', 'junction_support_pe', 'junction_support_sr', 'junction_support_su',
     'confidence_interval_pos_be1','confidence_interval_pos_95_be1', 'confidence_interval_end_be1','confidence_interval_end_95_be1',
     'confidence_interval_pos_be2','confidence_interval_pos_95_be2', 'confidence_interval_end_be2','confidence_interval_end_95_be2',
     'tum_support_pe','tum_support_sr', 'tum_support_su',
     'norm_support_pe','norm_support_sr', 'norm_support_su']]

#save lumpy events df if desired
#evdf.to_csv(outPath,sep='\t',header=True,index=False)