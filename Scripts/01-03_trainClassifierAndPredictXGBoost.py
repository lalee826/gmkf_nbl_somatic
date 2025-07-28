from xgboost import XGBClassifier
import pandas as pd
import numpy as np
import os
from sklearn.model_selection import KFold
from sklearn import preprocessing as preprocessing
from sklearn import metrics as metrics

workdir = '/rocker-build/gmkf_nbl_somatic/'

### Import and prepare training data
trainSet = workdir + 'Data/deepSVR_traindata.pkl'
train = pd.read_pickle(trainSet)
#remove AML samples
train = train[train['disease_AML'] != 1]
trainss = train.iloc[:10,:]
#check number of cases in each disease type
disTypes = ['disease_AML','disease_GST','disease_MPNST','disease_SCLC','disease_breast',
            'disease_colorectal','disease_glioblastoma','disease_lymphoma','disease_melanoma']
train.groupby(disTypes).size()
'''
counts of variants by disease types:
melanoma: 285
lymphoma: 1891
glioblastoma: 1256
crc: 1261
breast: 13306
SCLC: 13778
MPNST: 430
GST: 101
'''
#change germline label to fail
train.groupby(['call']).size()
''' #counts of different labels
a: ambiguous - 8864
f: fail - 5177
g: germline - 2516
s: somatic - 15751
'''
train['call'] = train['call'].replace('g','f')
train.groupby('call').size()
''' #new counts of different labels
a: ambiguous - 8864
f: fail - 7693
s: somatic - 15751
'''
#sort by columns
train.sort_index(axis=1, inplace=True)
train.columns.values
#take only the columns that we want (features that we could extract for our own data)
train = train.drop(['disease_AML', 'disease_GST', 'disease_MPNST',
       'disease_SCLC', 'disease_breast', 'disease_colorectal',
       'disease_glioblastoma', 'disease_lymphoma', 'disease_melanoma','reviewer_1', 'reviewer_2',
       'reviewer_3', 'reviewer_4'],axis=1)
train.shape # there are 58 features and 1 label
train = train.rename(columns={'normal_ref_avg_num_mismaches_as_fraction':'normal_ref_avg_num_mismatches_as_fraction',
                      'tumor_ref_avg_num_mismaches_as_fraction':'tumor_ref_avg_num_mismatches_as_fraction',
                      'normal_var_avg_num_mismaches_as_fraction':'normal_var_avg_num_mismatches_as_fraction',
                      'tumor_var_avg_num_mismaches_as_fraction':'tumor_var_avg_num_mismatches_as_fraction'})
trainCols = train.columns.values

### Split train data into features and labels
x_train = train.drop(['call'],axis=1).astype(float).values #get features
'''
0 - ambiguous
1 - fail
2 - somatic
'''
#min max scale train data
mm_scaler = preprocessing.MinMaxScaler()
x_train = mm_scaler.fit_transform(x_train)

#one-hot encode labels
y_train = pd.get_dummies(train.call).astype(int).values

#### CROSS VALIDATION OF BEST MODEL PARAMETERS ####
#store cross-val results
xgbCV = {}

#implement 5-fold split
kf = KFold(n_splits=5, shuffle=True)

#change y train array from one-hot encoded to single dimension
y_train_1d = np.argmax(y_train, axis=1)

for name,param in [('100,0.1,6,auc',{'n_estimators':100,'eta':0.1,'max_depth':6,'eval_metric':'auc'}),
                   ('100,0.3,10,auc',{'n_estimators':100,'eta':0.3,'max_depth':10,'eval_metric':'auc'}),
                   ('500,0.1,6,auc',{'n_estimators':500,'eta':0.1,'max_depth':6,'eval_metric':'auc'}),
                   ('500,0.3,10,auc',{'n_estimators':500,'eta':0.3,'max_depth':10,'eval_metric':'auc'}),
                   ('1000,0.1,6,auc',{'n_estimators':1000,'eta':0.1,'max_depth':6,'eval_metric':'auc'}),
                   ('1000,0.3,10,auc',{'n_estimators':1000,'eta':0.3,'max_depth':10,'eval_metric':'auc'}),
                   ('100,0.1,6,merror',{'n_estimators':100,'eta':0.1,'max_depth':6,'eval_metric':'merror'}),
                   ('100,0.3,10,merror',{'n_estimators':100,'eta':0.3,'max_depth':10,'eval_metric':'merror'}),
                   ('500,0.1,6,merror',{'n_estimators':500,'eta':0.1,'max_depth':6,'eval_metric':'merror'}),
                   ('500,0.3,10,merror',{'n_estimators':500,'eta':0.3,'max_depth':10,'eval_metric':'merror'}),
                   ('1000,0.1,6,merror',{'n_estimators':1000,'eta':0.1,'max_depth':6,'eval_metric':'merror'}),
                   ('1000,0.3,10,merror',{'n_estimators':1000,'eta':0.3,'max_depth':10,'eval_metric':'merror'})
                   ]:
    
    #initmodel
    xgbModel = XGBClassifier(**param)
    
    print(name)
    
    #store cv scores
    acc = []
    prec = {0:[],1:[],2:[]}
    recall = {0:[],1:[],2:[]}
    f1 = {0:[],1:[],2:[]}
    #spec = []
    
    #cross val
    for trainix,testix in kf.split(x_train):
        xtrain,ytrain = x_train[trainix,:], y_train_1d[trainix]
        xtest,ytest = x_train[testix,:], y_train_1d[testix]
        
        #fit model and predict
        xgbModel.fit(xtrain,ytrain)
        ypred = xgbModel.predict(xtest)        

        #evaluate performance
        from sklearn import metrics
        ac = round(metrics.accuracy_score(ytest,ypred),2)
        pr = metrics.precision_score(ytest,ypred,average=None).round(2)
        re = metrics.recall_score(ytest,ypred,average=None).round(2)
        f1s = metrics.f1_score(ytest,ypred,average=None).round(2)
        
        #save score for this kfold test
        acc.append(ac)
        prec[0].append(pr[0])
        prec[1].append(pr[1])
        prec[2].append(pr[2])
        recall[0].append(re[0])
        recall[1].append(re[1])
        recall[2].append(re[2])
        f1[0].append(f1s[0])
        f1[1].append(f1s[1])
        f1[2].append(f1s[2])
        
    #save parameter metrics
    xgbCV[name] = {'name':name,
     'parameters':param,
     'avgAccuracy':np.mean(acc).round(2),
     'avgPrecision_f':np.mean(prec[0]).round(2),
     'avgPrecision_a':np.mean(prec[1]).round(2),
     'avgPrecision_s':np.mean(prec[2]).round(2),
     'avgRecall_f':np.mean(recall[0]).round(2),
     'avgRecall_a':np.mean(recall[1]).round(2),
     'avgRecall_s':np.mean(recall[2]).round(2),
     'avgF1_f':np.mean(f1[0]).round(2),
     'avgF1_a':np.mean(f1[1]).round(2),
     'avgF1_s':np.mean(f1[2]).round(2)}   

###predict gmkf nbl data using best params
'''
Many parameter combinations performed equally well. We use the simplest combination.
'100,0.3,10,merror': {'avgAccuracy': 0.9,
  'avgF1_a': 0.93,
  'avgF1_f': 0.83,
  'avgF1_s': 0.92,
  'avgPrecision_a': 0.93,
  'avgPrecision_f': 0.86,
  'avgPrecision_s': 0.91,
  'avgRecall_a': 0.93,
  'avgRecall_f': 0.8,
  'avgRecall_s': 0.94,
  'name': '100,0.3,10,merror',
  'parameters': {'eta': 0.3,
   'eval_metric': 'merror',
   'max_depth': 10,
   'n_estimators': 100}
'''
#configure the model based on the data
xgbModel = XGBClassifier(n_estimators=100,eta=0.3,max_depth=10,eval_matric='merror')
xgbModel.fit(train.drop('call',axis=1),y_train_1d)

### import our own data
consensusVarsFile = workdir + 'Data/coding_mutations_annotated.tsv'
uniqueVarsFile = workdir + 'Data/unique_coding_mutations_full.tsv'
bamReadCountsFiles = workdir + 'Data/variants_bamrc_analyzed.tsv'

consensusVars = pd.read_table(consensusVarsFile,encoding="latin1")
uniqueVars = pd.read_table(uniqueVarsFile)
brcData = pd.read_table(bamReadCountsFiles)

#merge consensus variants and unique variants and grab bam-readcount output
nblData = pd.concat([brcData.merge(consensusVars[['mutid']],on='mutid',how='inner'),
                   brcData.merge(uniqueVars[['mutid']],on='mutid',how='inner')]).reset_index(drop=True)
#set mutid as index
nblData.set_index('mutid',inplace=True)
#rename columns to match train data
nblData = nblData.rename({'normal_ref_avg_num_mismaches_as_fraction':'normal_ref_avg_num_mismatches_as_fraction',
              'tumor_ref_avg_num_mismaches_as_fraction':'tumor_ref_avg_num_mismatches_as_fraction',
              'normal_var_avg_num_mismaches_as_fraction':'normal_var_avg_num_mismatches_as_fraction',
              'tumor_var_avg_num_mismaches_as_fraction':'tumor_var_avg_num_mismatches_as_fraction'},axis=1)

#order columns to match train data
matchCols = np.delete(trainCols,np.where(trainCols=='call')[0])
nblData = nblData[matchCols]

#make predictions using trained model
predLabelsNBL = xgbModel.predict(nblData)

#look at prediction counts
print(np.bincount(predLabelsNBL))

#add classifier predictions to each variant
nblData['prediction'] = predLabelsNBL
nblData['prediction'].replace({0:'ambiguous',1:'fail',2:'somatic'},inplace=True)

#join back to consensus and unique variant mutations files
consensusVars = consensusVars.merge(nblData,on='mutid',how='left')
uniqueVars = uniqueVars.merge(nblData,on='mutid',how='left')