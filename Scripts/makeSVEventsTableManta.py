"""
Description: Parse SV calls from Manta to create a table of SV events. This script will:
    
    - Join Breakpoints into events with multiple ends
    - Tabulate significant metrics
    - Provide a tabular view of final somatic filtered SV calls made by Manta
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
df = pd.read_table('/rocker-build/gmkf_nbl_somatic/Raw Data/Manta VCF/PATBFN.manta.clean.vcf',sep='\t',header=0)

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
#csq keys
csqKeys = ('Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|'
      'CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|PICK|'
      'VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|RefSeq|GENE_PHENO|'
      'SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|'
      'gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|'
      'gnomAD_SAS_AF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE').split('|')

###extract data from info and tumor/normal columns
info = df.INFO
formats = df.FORMAT
tumorValues = df.iloc[:,tumor_col]
normalValues = df.iloc[:,normal_col]

###create these to make table of sv calls
info_dicts = [] #metrics for variant call in list of dictionaries
csq_dicts = [] #list of csq metrics for each call
is_imprecise = [] #bool for IMPRECISE metric
is_somatic = [] #bool for SECONDARY metric
tnv = [] #metrics for tumor and normal in list of dictionaries


for i in np.arange(0,len(info)):
    
    #i = 9
    
    #split info field
    infoSplit = info[i].split(';')
    #split each info field by key,value
    infoSplit = [x.split('=') for x in infoSplit]
    #check if call is imprecise
    isImprecise = np.any([x[0]=='IMPRECISE' for x in infoSplit])
    #check if call is secondary
    isSomatic = np.any([x[0]=='SOMATIC' for x in infoSplit])
    #get csq fields and remove from info
    csqIX = np.arange(0,len(infoSplit))[np.array([x[0] == 'CSQ' for x in infoSplit])][0]
    csq = infoSplit[csqIX]
    del infoSplit[csqIX]
    #ensure each csq entry is of same length as keys
    csqVals = csq[1].split(',')
    for entry in csqVals:
        try:
            assert len(entry.split('|')) == len(csqKeys)
        except: 
            raise AssertionError(str(i)+':'+entry)
    #save csq values by key and value to list
    csq_dicts_sample = [dict(zip(csqKeys,csqVals[i].split('|'))) for i in np.arange(0,len(csqVals))]
    #save csq list
    csq_dicts.append(csq_dicts_sample)
    
    '''
    #remove imprecise and check if any other fields are len==1
    if np.any([len(y) != 2 for y in list(np.array(infoSplit)[[x[0] not in ['IMPRECISE','SECONDARY'] for x in infoSplit]])]):
        print(i)
    '''
    
    #remove any fields that are not of key=value format
    infoSplit = list(np.array(infoSplit,dtype=object)[[len(x)==2 for x in infoSplit]])
    #create dictionary of fields and values
    #save row to lists
    info_dicts.append(dict(tuple([tuple(x) for x in infoSplit])))
    is_imprecise.append(isImprecise)
    is_somatic.append(isSomatic)
    
    #get field keys
    form = formats[i].split(':')
    #get values for tumor and normal
    tv = tumorValues[i].split(':')
    nv = normalValues[i].split(':')
    #create tumor and normal dictionaries
    td = dict(tuple([tuple([form[j],tv[j]]) for j in np.arange(0,len(form))]))
    nd = dict(tuple([tuple([form[j],nv[j]]) for j in np.arange(0,len(form))]))
    tnv.append({'tum':td,'norm':nd})

#%% Extract the values to be included in table
######### PRIMARY METRICS
#pull chroms
chroms = df.CHROM.tolist()
#pull POS
pos = df.POS.tolist()
#pull end position
ends = [info_dicts[i]['END'] if 'END' in info_dicts[i] else 'NaN' for i in np.arange(0,len(info_dicts))]
#pull sv length
svlen = [abs(int(info_dicts[i]['SVLEN'])) if 'SVLEN' in info_dicts[i] else 'NaN' for i in np.arange(0,len(info_dicts))]
#pull insertion length and sequence
svinsseq = [info_dicts[i]['SVINSSEQ'] if 'SVINSSEQ' in info_dicts[i] else 'NaN' for i in np.arange(0,len(info_dicts))]
svinslen = [abs(int(info_dicts[i]['SVINSLEN'])) if 'SVINSLEN' in info_dicts[i] else 'NaN' for i in np.arange(0,len(info_dicts))]
#pull ID
ids = df.ID.tolist() #manta call IDs
#pull event
events = [x[:-1] for x in ids]
#get snv type
types = [x.split(':')[0][-3:] for x in ids]
#get filter metric
fil = df.FILTER.tolist()

######### TUMOR NORMAL BREAKEND METRICS
##get spanning split and paired end reads for ref/alt in tumor and normal
t = [tnv[i]['tum']['SR'].split(',')  if 'SR' in tnv[i]['tum'] else 'NaN' for i in np.arange(0,len(info_dicts))]
trefsr = [x if x == 'NaN' else x[0] for x in t]
taltsr = [x if x == 'NaN' else x[1] for x in t]
t = [tnv[i]['tum']['PR'].split(',')  if 'PR' in tnv[i]['tum'] else 'NaN' for i in np.arange(0,len(info_dicts))]
trefpr = [x if x == 'NaN' else x[0] for x in t]
taltpr = [x if x == 'NaN' else x[1] for x in t]

t = [tnv[i]['norm']['SR'].split(',')  if 'SR' in tnv[i]['norm'] else 'NaN' for i in np.arange(0,len(info_dicts))]
nrefsr = [x if x == 'NaN' else x[0] for x in t]
naltsr = [x if x == 'NaN' else x[1] for x in t]
t = [tnv[i]['norm']['PR'].split(',')  if 'PR' in tnv[i]['norm'] else 'NaN' for i in np.arange(0,len(info_dicts))]
nrefpr = [x if x == 'NaN' else x[0] for x in t]
naltpr = [x if x == 'NaN' else x[1] for x in t]

###################pull from info tags
#get all possible info tags
t = [list(info_dicts[i].keys()) for i in np.arange(0,len(info_dicts))]
info_tags = set([item for sublist in t for item in sublist])
bndd = [info_dicts[i]['BND_DEPTH'] if 'BND_DEPTH' in info_dicts[i] else 'NaN' for i in np.arange(0,len(info_dicts))]
ciend = [info_dicts[i]['CIEND'] if 'CIEND' in info_dicts[i] else 'NaN' for i in np.arange(0,len(info_dicts))]
cigar = [info_dicts[i]['CIGAR'] if 'CIGAR' in info_dicts[i] else 'NaN' for i in np.arange(0,len(info_dicts))]
homlen = [info_dicts[i]['HOMLEN'] if 'HOMLEN' in info_dicts[i] else 'NaN' for i in np.arange(0,len(info_dicts))]
homseq = [info_dicts[i]['HOMSEQ'] if 'HOMSEQ' in info_dicts[i] else 'NaN' for i in np.arange(0,len(info_dicts))]
mid = [info_dicts[i]['MATEID'] if 'MATEID' in info_dicts[i] else 'NaN' for i in np.arange(0,len(info_dicts))]
somscore = [info_dicts[i]['SOMATICSCORE'] if 'SOMATICSCORE' in info_dicts[i] else 'NaN' for i in np.arange(0,len(info_dicts))]

######### CSQ METRICS
'''   Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|
      CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|PICK|
      VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|RefSeq|GENE_PHENO|
      SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|
      gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|
      gnomAD_SAS_AF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE'''
csq_allele = [('|').join(k) for k in [[j ['Allele'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_consequence = [('|').join(k) for k in [[j ['Consequence'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_impact = [('|').join(k) for k in [[j ['IMPACT'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_symbol = [('|').join(k) for k in [[j ['SYMBOL'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_gene = [('|').join(k) for k in [[j ['Gene'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_featuretype = [('|').join(k) for k in [[j ['Feature_type'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_feature = [('|').join(k) for k in [[j ['Feature'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_biotype = [('|').join(k) for k in [[j ['BIOTYPE'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_exon = [('|').join(k) for k in [[j ['EXON'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_intron = [('|').join(k) for k in [[j ['INTRON'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_hgvsc = [('|').join(k) for k in [[j ['HGVSc'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_cdnapos = [('|').join(k) for k in [[j ['cDNA_position'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_cdspos = [('|').join(k) for k in [[j ['CDS_position'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_proteinpos = [('|').join(k) for k in [[j ['Protein_position'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_distance = [('|').join(k) for k in [[j ['DISTANCE'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_strand = [('|').join(k) for k in [[j ['STRAND'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_flags = [('|').join(k) for k in [[j ['FLAGS'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_pick = [('|').join(k) for k in [[j ['PICK'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_symbolsource = [('|').join(k) for k in [[j ['SYMBOL_SOURCE'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_hgncid = [('|').join(k) for k in [[j ['HGNC_ID'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_canonical = [('|').join(k) for k in [[j ['CANONICAL'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_ccds = [('|').join(k) for k in [[j ['CCDS'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_ensp = [('|').join(k) for k in [[j ['ENSP'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_swissprot = [('|').join(k) for k in [[j ['SWISSPROT'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_trembl = [('|').join(k) for k in [[j ['TREMBL'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_uniparc = [('|').join(k) for k in [[j ['UNIPARC'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_refseq = [('|').join(k) for k in [[j ['RefSeq'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_domains = [('|').join(k) for k in [[j ['DOMAINS'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_af = [('|').join(k) for k in [[j ['AF'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_gnomadaf = [('|').join(k) for k in [[j ['gnomAD_AF'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_clinsig = [('|').join(k) for k in [[j ['CLIN_SIG'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]
csq_somatic = [('|').join(k) for k in [[j ['SOMATIC'] for j in csq_dicts[i]] for i in np.arange(0,len(csq_dicts))]]

#%% create a data frame to output
svdf = pd.DataFrame({'event':events,'id':ids,'mateid':mid,'type':types,'chrom':chroms,'pos':pos,'end':ends,'length':svlen, #primary metrics
                     'imprecise':is_imprecise,'somatic':is_somatic,'filter':fil, #primary metrics
                     'insertion_length':svinslen, 'insertion_sequence':svinsseq, #insertion metrics
                     'bnd_depth':bndd,'confidence_interval_end':ciend,'cigar':cigar,'homology_length':homlen,#info tags
                     'homology_sequence':homseq,'somatic_score':somscore,#info tags
                     'tumor_alt_pe':taltpr,'tumor_alt_sr':taltsr,'tumor_ref_pe':trefpr,'tumor_ref_sr':trefsr, #tumor support metrics
                     'normal_alt_pe':naltpr,'normal_alt_sr':naltsr,'normal_ref_pe':nrefpr,'normal_ref_sr':nrefsr, #normal support metrics
                     'csq_allele':csq_allele,'csq_consequence':csq_consequence,'csq_impact':csq_impact,'csq_symbol':csq_symbol, #csq metrics
                     'csq_gene':csq_gene,'csq_featuretype':csq_featuretype,'csq_feature':csq_feature,'csq_biotype':csq_biotype,
                     'csq_exon':csq_exon,'csq_intron':csq_intron,'csq_hgvsc':csq_hgvsc,'csq_cdnapos':csq_cdnapos,
                     'csq_cdspos':csq_cdspos,'csq_proteinpos':csq_proteinpos,'csq_distance':csq_distance,'csq_strand':csq_strand,
                     'csq_flags':csq_flags,'csq_pick':csq_pick,'csq_symbolsource':csq_symbolsource,'csq_hgncid':csq_hgncid,
                     'csq_canonical':csq_canonical,'csq_ccds':csq_ccds,'csq_ensp':csq_ensp,'csq_swissprot':csq_swissprot,
                     'csq_trembl':csq_trembl,'csq_uniparc':csq_uniparc,'csq_refseq':csq_refseq,'csq_domains':csq_domains,
                     'csq_af':csq_af,'csq_gnomadaf':csq_gnomadaf,'csq_clinsig':csq_clinsig,'csq_somatic':csq_somatic})
#sort by chromosome
svdf.chrom = pd.Categorical(svdf.chrom, 
                      categories=['chr'+str(x) for x in list(range(1,23))+['X','Y']],
                      ordered=True)

#%% combine mate pairs and make events df
#### For BND events join mate pairs, discard those without mates
#split off BND from rest
svdf_bnd = svdf.loc[svdf.type=='BND',:].reset_index(drop=True)
svdf = svdf.loc[svdf.type!='BND',:].reset_index(drop=True)
#drop mate id from non bnd structural variants
svdf = svdf.drop(['mateid'],axis=1)
#discard rows whose mate pair is not in data frame
allID = svdf_bnd.id.tolist()
svdf_bnd = svdf_bnd.loc[[i in allID for i in svdf_bnd.mateid],:].reset_index(drop=True)

#make list of columns that are kept separate between two breakpoints in each event
split_cols = ['id','chrom','pos','imprecise','somatic','filter','bnd_depth','homology_length','homology_sequence',
              'csq_allele', 'csq_consequence','csq_impact','csq_symbol','csq_gene','csq_featuretype',
              'csq_feature','csq_biotype','csq_exon','csq_intron','csq_hgvsc','csq_cdnapos','csq_cdspos',
              'csq_proteinpos','csq_distance','csq_strand','csq_flags','csq_pick','csq_symbolsource',
              'csq_hgncid','csq_canonical','csq_ccds','csq_ensp','csq_swissprot','csq_trembl','csq_uniparc',
              'csq_refseq','csq_domains','csq_af','csq_gnomadaf','csq_clinsig','csq_somatic']
common_cols_events = ['event','type','length','insertion_length','insertion_sequence',
                      'confidence_interval_end','cigar','somatic_score','tumor_alt_pe','tumor_alt_sr',
                      'tumor_ref_pe','tumor_ref_sr','normal_alt_pe','normal_alt_sr','normal_ref_pe','normal_ref_sr']

#make new column names for values that are kept separate when breakends are combined into an event
split_cols_events = [[x+sffx for sffx in ['_be1','_be2']] for x in split_cols]
split_cols_events = [y for x in split_cols_events for y in x]
#make list of column names for events df
edf_cols = common_cols_events+split_cols_events

#convert to events df
evdf = pd.DataFrame()
for ix,ev in enumerate(set(svdf_bnd.event)):
    '''
    ix = 0
    ev = 'MantaBND:143037:0:1:0:0:0:'
    '''
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

#%% combine the events df with non bnd svs
#modify non bnd columns and names and bind with BND events
rename_dict = dict(zip(split_cols,[x+'_be1' for x in split_cols]))
rename_dict['end'] = 'pos_be2'
svdf = svdf.rename(columns=rename_dict)
#add NaN values for be2 columns
nanFill = ['NaN']*svdf.shape[0]
new_cols_dict = dict(zip([x+'_be2' for x in split_cols],[nanFill]*len(split_cols)))
#don't change pos_be2
del new_cols_dict['pos_be2']
#chromosome on be2 is same as be1
new_cols_dict['chrom_be2'] = svdf.chrom_be1.tolist()
#add new columns to non bnd structural variants
svdf = pd.concat([svdf,pd.DataFrame(new_cols_dict)],axis=1)

###append bnd snvs to rest
evdf = svdf.append(evdf,ignore_index=True)
#order by breakend 1 chrom and pos
evdf = evdf.sort_values(['chrom_be1','pos_be1']).reset_index(drop=True)
###reorder columns
#evdf.columns.values
col_order = ['event','id_be1','id_be2','type','chrom_be1','pos_be1','chrom_be2','pos_be2','length','imprecise_be1','imprecise_be2',
     'somatic_score','somatic_be1','somatic_be2','filter_be1','filter_be2','insertion_length','insertion_sequence',
     'bnd_depth_be1','bnd_depth_be2','confidence_interval_end','cigar','homology_length_be1','homology_length_be2',
     'homology_sequence_be1','homology_sequence_be2',
     'tumor_alt_pe','tumor_alt_sr','tumor_ref_pe','tumor_ref_sr','normal_alt_pe','normal_alt_sr','normal_ref_pe','normal_ref_sr']
evdf = evdf[col_order]
### sort and change data types to integers
#sort by chromosome
evdf.chrom_be1 = pd.Categorical(evdf.chrom_be1, 
                      categories=['chr'+str(x) for x in list(range(1,23))+['X','Y']],
                      ordered=True)
evdf.chrom_be2 = pd.Categorical(evdf.chrom_be2, 
                      categories=['chr'+str(x) for x in list(range(1,23))+['X','Y']],
                      ordered=True)
evdf = evdf.sort_values(['chrom_be1','pos_be1','chrom_be2','pos_be2']).reset_index(drop=True)
evdf.pos_be1 = evdf.pos_be1.astype(int)
evdf.pos_be2 = evdf.pos_be2.astype(int)
evdf.imprecise_be1 = evdf.imprecise_be1.astype(bool)
evdf.imprecise_be2 = evdf.imprecise_be2.astype(bool)
evdf.somatic_be1 = evdf.somatic_be1.astype(bool)
evdf.somatic_be2 = evdf.somatic_be2.astype(bool)

#Save Manta events table if desired
#evdf.to_csv(outPath,sep='\t',header=True,index=False)