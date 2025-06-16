'''
- This script evaluates different classifier models on DeepSVR training data tailored to our data set
- We remove liquid tumor variants and collapse the feature set to 50 extractable features from our own data set
- Two data files are necessary for this script:
    1. DeepSVR training data
    2. Bam-Readcount features extracted for all variants in our data (provided in repository without script)
'''

#%%Import packages
import pandas as pd
import numpy as np
import os
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from sklearn import metrics
import matplotlib.pyplot as plt
from itertools import cycle

workdir = '/rocker-build/gmkf_nbl_somatic/'

#%% Import and prepare training data
#Import DeepSVR training data and bam-readcount features extracted for our variants
trainFile = workdir + 'Data/DeepSVR_training_data_preprocessed.pkl'
dataFile = workdir + 'Data/variants_bamrc_analyzed.tsv'
data = pd.read_table(dataFile)
train = pd.read_pickle(trainFile)
train.shape #(41832, 72)
#remove AML samples
train = train[train['disease_AML'] != 1]
train.shape #(32308, 72)
train['call'] = train['call'].replace('g','f') #germline samples should be considered fail
#sort by columns
train.sort_index(axis=1, inplace=True)
train.columns.values
train.shape #(32308, 72)
#take only the columns that we want (features that we could extract for our own data)
train = train.drop(['disease_AML', 'disease_GST', 'disease_MPNST',
       'disease_SCLC', 'disease_breast', 'disease_colorectal',
       'disease_glioblastoma', 'disease_lymphoma', 'disease_melanoma','reviewer_1', 'reviewer_2',
       'reviewer_3', 'reviewer_4'],axis=1)
train.shape #(32308, 59)

#index our data by mutid and order
data.index = data['mutid'].tolist()
var_order = data.index.values
data = data.drop(['mutid'],axis=1)

#rename some typos in column names
data = data.rename(columns={'normal_ref_avg_num_mismaches_as_fraction':'normal_ref_avg_num_mismatches_as_fraction',
                            'tumor_var_avg_num_mismaches_as_fraction':'tumor_var_avg_num_mismatches_as_fraction'})
train = train.rename(columns={'normal_ref_avg_num_mismaches_as_fraction':'normal_ref_avg_num_mismatches_as_fraction',
                              'normal_var_avg_num_mismaches_as_fraction':'normal_var_avg_num_mismatches_as_fraction',
                              'tumor_ref_avg_num_mismaches_as_fraction':'tumor_ref_avg_num_mismatches_as_fraction',
                              'tumor_var_avg_num_mismaches_as_fraction':'tumor_var_avg_num_mismatches_as_fraction'})
col_order = data.columns.values.tolist() + ['call']
#make column order consistent between training data and 
train = train[col_order]
train.columns.values
#min-max scale columns in train data
mm_scaler = preprocessing.MinMaxScaler()
train_scaled = mm_scaler.fit_transform(train.drop(['call'],axis=1))
train_unscaled = train.drop(['call'],axis=1)

#%% Split into train/test sets, and then cross-validation sets
Y = train.call.replace({'a':0,'f':1,'s':2}).astype(int)
#X = train.drop(['call'],axis=1).astype(float).values #get features
X = train_unscaled
#fix random seed for reproducibility
seed = 999
np.random.seed(seed)
#split into train/test
x_train,x_test,y_train,y_test = train_test_split(X,Y,test_size=0.1,random_state=seed)
#set cross-validation params
kfold = KFold(n_splits=5,shuffle=True,random_state=seed)
### test best params for random forest model in cross-validation
fold_acc = {}
fold_prec = {}
fold_rec = {}
fold_f1 = {}
cores = 4 #parallel jobs

for n in [10,50,100,500,1000]:

    print(n)
    
    cv_acc = []
    cv_prec = {'a':[],'f':[],'s':[]}
    cv_rec = {'a':[],'f':[],'s':[]}
    cv_f1 = {'a':[],'f':[],'s':[]}

    #init model
    rfc = RandomForestClassifier(n_estimators=n,criterion='gini',max_depth=None,min_samples_split=2,min_samples_leaf=1,
                             min_weight_fraction_leaf=0.0,max_features='auto',max_leaf_nodes=None,min_impurity_decrease=0.0,min_impurity_split=None,
                             bootstrap=True,oob_score=False,n_jobs=cores,random_state=None,verbose=0,warm_start=False,class_weight=None)
    
    for trainix,testix in kfold.split(X=train_unscaled,y=Y):
        
        #get train,test sets
        x_train = X.iloc[trainix]
        y_train = Y.iloc[trainix]
        x_test = X.iloc[testix]
        y_test = Y.iloc[testix]
        
        #fit model and predict
        rfc.fit(x_train,y_train)
        pred_labels = rfc.predict(x_test)
    
        #assess performance
        acc = round(metrics.accuracy_score(y_test, pred_labels),2)
        ps = metrics.precision_score(y_test,pred_labels,average=None)
        rs = metrics.recall_score(y_test,pred_labels,average=None)
        fs = metrics.f1_score(y_test,pred_labels,average=None)
    
        #save values
        cv_acc.append(acc)
        cv_prec['a'] = cv_prec['a'] + [round(ps[0],2)]
        cv_prec['f'] = cv_prec['f'] + [round(ps[1],2)]
        cv_prec['s'] = cv_prec['s'] + [round(ps[2],2)]
        cv_rec['a'] = cv_rec['a'] + [round(rs[0],2)]
        cv_rec['f'] = cv_rec['f'] + [round(rs[1],2)]
        cv_rec['s'] = cv_rec['s'] + [round(rs[2],2)]
        cv_f1['a'] = cv_f1['a'] + [round(fs[0],2)]
        cv_f1['f'] = cv_f1['f'] + [round(fs[1],2)]
        cv_f1['s'] = cv_f1['s'] + [round(fs[2],2)]
    
    #get average cross-train scores
    avg_acc = round(np.mean(cv_acc),2)
    avg_prec = {'a':round(np.mean(cv_prec['a']),2),'f':round(np.mean(cv_prec['f']),2),'s':round(np.mean(cv_prec['s']),2)}
    avg_rec = {'a':round(np.mean(cv_rec['a']),2),'f':round(np.mean(cv_rec['f']),2),'s':round(np.mean(cv_rec['s']),2)}
    avg_f1 = {'a':round(np.mean(cv_f1['a']),2),'f':round(np.mean(cv_f1['f']),2),'s':round(np.mean(cv_f1['s']),2)}
    #save cross-train scores for iteration
    fold_acc[n] = avg_acc
    fold_prec[n] = avg_prec
    fold_rec[n] = avg_rec  
    fold_f1[n] = avg_f1    

#%%test best params for RBM kernel SVM model
#scale the data using standard scaler
std_scaler = preprocessing.StandardScaler()
train_scaled = std_scaler.fit_transform(train.drop(['call'],axis=1))
### Split into train/test sets, and then cross-validation sets
Y = train.call.replace({'a':0,'f':1,'s':2}).astype(float)
#X = train.drop(['call'],axis=1).astype(float).values #get features
X = train_scaled.astype(float)
#fix random seed for reproducibility
seed = 999
np.random.seed(seed)
#split into train/test
x_train,x_test,y_train,y_test = train_test_split(X,Y,test_size=0.1,random_state=seed)
#set cross-validation params
kfold = KFold(n_splits=5,shuffle=True,random_state=seed)
###test best regularization constant
fold_acc = {}
fold_prec = {}
fold_rec = {}
fold_f1 = {}

for c in [0.000001,0.00001,0.0001,0.001,0.01,0.1,1,10,100]:
    
    print(c)
    
    cv_acc = []
    cv_prec = {'a':[],'f':[],'s':[]}
    cv_rec = {'a':[],'f':[],'s':[]}
    cv_f1 = {'a':[],'f':[],'s':[]}
    
    #init model
    svmc = svm.SVC(C=c, kernel='rbf',gamma='auto',shrinking=True,probability=False,tol=0.001,cache_size=200,class_weight=None,
                   verbose=0,max_iter=-1,random_state=None)
    
    #cross-validation
    for trainix,testix in kfold.split(X=X,y=Y):
        
        #split into train,test sets
        x_train = X[trainix]
        y_train = Y.iloc[trainix]
        x_test = X[testix]
        y_test = Y.iloc[testix]
        
        #fit model and predict
        svmc.fit(x_train,y_train)
        pred_labels = svmc.predict(x_test)
    
        #assess performance
        acc = round(metrics.accuracy_score(y_test, pred_labels),2)
        ps = metrics.precision_score(y_test,pred_labels,average=None)
        rs = metrics.recall_score(y_test,pred_labels,average=None)
        fs = metrics.f1_score(y_test,pred_labels,average=None)
    
        #save values
        cv_acc.append(acc)
        cv_prec['a'] = cv_prec['a'] + [round(ps[0],2)]
        cv_prec['f'] = cv_prec['f'] + [round(ps[1],2)]
        cv_prec['s'] = cv_prec['s'] + [round(ps[2],2)]
        cv_rec['a'] = cv_rec['a'] + [round(rs[0],2)]
        cv_rec['f'] = cv_rec['f'] + [round(rs[1],2)]
        cv_rec['s'] = cv_rec['s'] + [round(rs[2],2)]
        cv_f1['a'] = cv_f1['a'] + [round(fs[0],2)]
        cv_f1['f'] = cv_f1['f'] + [round(fs[1],2)]
        cv_f1['s'] = cv_f1['s'] + [round(fs[2],2)]
    
    #get average cross-train scores
    avg_acc = round(np.mean(cv_acc),2)
    avg_prec = {'a':round(np.mean(cv_prec['a']),2),'f':round(np.mean(cv_prec['f']),2),'s':round(np.mean(cv_prec['s']),2)}
    avg_rec = {'a':round(np.mean(cv_rec['a']),2),'f':round(np.mean(cv_rec['f']),2),'s':round(np.mean(cv_rec['s']),2)}
    avg_f1 = {'a':round(np.mean(cv_f1['a']),2),'f':round(np.mean(cv_f1['f']),2),'s':round(np.mean(cv_f1['s']),2)}
    #save cross-train scores for iteration
    fold_acc[c] = avg_acc
    fold_prec[c] = avg_prec
    fold_rec[c] = avg_rec  
    fold_f1[c] = avg_f1    
    
    print(avg_acc)

#%% find best params for gradient boosting classifier
Y = train.call.replace({'a':0,'f':1,'s':2}).astype(int)
#X = train.drop(['call'],axis=1).astype(float).values #get features
X = train_unscaled
#fix random seed for reproducibility
seed = 999
np.random.seed(seed)
#split into train/test
x_train,x_test,y_train,y_test = train_test_split(X,Y,test_size=0.1,random_state=seed)
#set cross-validation params
kfold = KFold(n_splits=5,shuffle=True,random_state=seed)
### test best params for random forest model in cross-validation

''' default settings
l = 0.1
n = 100
d = 3
'''

params = {'learning_rate':0.1,'n_estimators':100,'max_depth':3,'loss':'deviance','subsample':1.0,'min_samples_split':2,
                                 'min_samples_leaf':1,'min_weight_fraction_leaf':0.0,'min_impurity_decrease':0.0,
                                 'verbose':0,'max_leaf_nodes':None}
#test parameters in k-fold cross-validation
param_acc = {}
param_prec = {}
param_rec = {}
param_f1 = {}
for name,setting in [('0.001,100,3',{'learning_rate': 0.001,'n_estimators':100,'max_depth':3}),
                     ('0.01,100,3',{'learning_rate': 0.01,'n_estimators':100,'max_depth':3}),
                     ('0.1,100,3',{'learning_rate': 0.1,'n_estimators':100,'max_depth':3}),
                     ('1,100,3',{'learning_rate': 1,'n_estimators':100,'max_depth':3}),
                     ('0.1,10,3',{'learning_rate': 0.1,'n_estimators':10,'max_depth':3}),
                     ('0.1,100,3',{'learning_rate': 0.1,'n_estimators':100,'max_depth':3}),
                     ('0.1,500,3',{'learning_rate': 0.1,'n_estimators':500,'max_depth':3}),
                     ('0.1,1000,3',{'learning_rate': 0.1,'n_estimators':1000,'max_depth':3}),
                     ('0.1,100,2',{'learning_rate': 0.1,'n_estimators':100,'max_depth':2}),
                     ('0.1,100,3',{'learning_rate': 0.1,'n_estimators':100,'max_depth':3}),
                     ('0.1,100,4',{'learning_rate': 0.1,'n_estimators':100,'max_depth':4}),
                     ('0.1,100,5',{'learning_rate': 0.1,'n_estimators':100,'max_depth':5})]:

    params.update(setting)
    print(name,params)

    #init model
    gbc = GradientBoostingClassifier(**params)

    #cross-validation
    cv_acc = []
    cv_prec = {'a':[],'f':[],'s':[]}
    cv_rec = {'a':[],'f':[],'s':[]}
    cv_f1 = {'a':[],'f':[],'s':[]}
    for trainix,testix in kfold.split(X=X,y=Y):
        
        #split into train,test sets
        x_train = X.iloc[trainix]
        y_train = Y.iloc[trainix]
        x_test = X.iloc[testix]
        y_test = Y.iloc[testix]
        
        #fit model and predict
        gbc.fit(x_train,y_train)
        pred_labels = gbc.predict(x_test)
        
        #assess performance
        acc = round(metrics.accuracy_score(y_test, pred_labels),2)
        ps = metrics.precision_score(y_test,pred_labels,average=None)
        rs = metrics.recall_score(y_test,pred_labels,average=None)
        fs = metrics.f1_score(y_test,pred_labels,average=None)
    
        #save values
        cv_acc.append(acc)
        cv_prec['a'] = cv_prec['a'] + [round(ps[0],2)]
        cv_prec['f'] = cv_prec['f'] + [round(ps[1],2)]
        cv_prec['s'] = cv_prec['s'] + [round(ps[2],2)]
        cv_rec['a'] = cv_rec['a'] + [round(rs[0],2)]
        cv_rec['f'] = cv_rec['f'] + [round(rs[1],2)]
        cv_rec['s'] = cv_rec['s'] + [round(rs[2],2)]
        cv_f1['a'] = cv_f1['a'] + [round(fs[0],2)]
        cv_f1['f'] = cv_f1['f'] + [round(fs[1],2)]
        cv_f1['s'] = cv_f1['s'] + [round(fs[2],2)]
        
    #get average cross-train scores
    avg_acc = round(np.mean(cv_acc),2)
    avg_prec = {'a':round(np.mean(cv_prec['a']),2),'f':round(np.mean(cv_prec['f']),2),'s':round(np.mean(cv_prec['s']),2)}
    avg_rec = {'a':round(np.mean(cv_rec['a']),2),'f':round(np.mean(cv_rec['f']),2),'s':round(np.mean(cv_rec['s']),2)}
    avg_f1 = {'a':round(np.mean(cv_f1['a']),2),'f':round(np.mean(cv_f1['f']),2),'s':round(np.mean(cv_f1['s']),2)}
    #save cross-train scores for iteration
    param_acc[name] = avg_acc
    param_prec[name] = avg_prec
    param_rec[name] = avg_rec  
    param_f1[name] = avg_f1
#### 2nd round of testing
#test parameters in k-fold cross-validation
param_acc_2 = {}
param_prec_2 = {}
param_rec_2 = {}
param_f1_2 = {}

'''
check 250,500,1000 estimators vs. 4-6 max depth
'''
for name,setting in [('0.1,250,4',{'learning_rate': 0.1,'n_estimators':250,'max_depth':4}),
                     ('0.1,500,4',{'learning_rate': 0.1,'n_estimators':500,'max_depth':4}),
                     ('0.1,1000,4',{'learning_rate': 0.1,'n_estimators':1000,'max_depth':4}),
                     ('0.1,250,5',{'learning_rate': 0.1,'n_estimators':250,'max_depth':5}),
                     ('0.1,500,5',{'learning_rate': 0.1,'n_estimators':500,'max_depth':5}),
                     ('0.1,1000,5',{'learning_rate': 0.1,'n_estimators':1000,'max_depth':5}),
                     ('0.1,250,6',{'learning_rate': 0.1,'n_estimators':250,'max_depth':6}),
                     ('0.1,500,6',{'learning_rate': 0.1,'n_estimators':500,'max_depth':6}),
                     ('0.1,1000,6',{'learning_rate': 0.1,'n_estimators':1000,'max_depth':6})
                     ]:

    params.update(setting)
    print(name,params)

    #init model
    gbc = GradientBoostingClassifier(**params)

    #cross-validation
    cv_acc = []
    cv_prec = {'a':[],'f':[],'s':[]}
    cv_rec = {'a':[],'f':[],'s':[]}
    cv_f1 = {'a':[],'f':[],'s':[]}
    for trainix,testix in kfold.split(X=X,y=Y):
        
        #split into train,test sets
        x_train = X.iloc[trainix]
        y_train = Y.iloc[trainix]
        x_test = X.iloc[testix]
        y_test = Y.iloc[testix]
        
        #fit model and predict
        gbc.fit(x_train,y_train)
        pred_labels = gbc.predict(x_test)
        
        #assess performance
        acc = round(metrics.accuracy_score(y_test, pred_labels),2)
        ps = metrics.precision_score(y_test,pred_labels,average=None)
        rs = metrics.recall_score(y_test,pred_labels,average=None)
        fs = metrics.f1_score(y_test,pred_labels,average=None)
    
        #save values
        cv_acc.append(acc)
        cv_prec['a'] = cv_prec['a'] + [round(ps[0],2)]
        cv_prec['f'] = cv_prec['f'] + [round(ps[1],2)]
        cv_prec['s'] = cv_prec['s'] + [round(ps[2],2)]
        cv_rec['a'] = cv_rec['a'] + [round(rs[0],2)]
        cv_rec['f'] = cv_rec['f'] + [round(rs[1],2)]
        cv_rec['s'] = cv_rec['s'] + [round(rs[2],2)]
        cv_f1['a'] = cv_f1['a'] + [round(fs[0],2)]
        cv_f1['f'] = cv_f1['f'] + [round(fs[1],2)]
        cv_f1['s'] = cv_f1['s'] + [round(fs[2],2)]
        
    #get average cross-train scores
    avg_acc = round(np.mean(cv_acc),2)
    avg_prec = {'a':round(np.mean(cv_prec['a']),2),'f':round(np.mean(cv_prec['f']),2),'s':round(np.mean(cv_prec['s']),2)}
    avg_rec = {'a':round(np.mean(cv_rec['a']),2),'f':round(np.mean(cv_rec['f']),2),'s':round(np.mean(cv_rec['s']),2)}
    avg_f1 = {'a':round(np.mean(cv_f1['a']),2),'f':round(np.mean(cv_f1['f']),2),'s':round(np.mean(cv_f1['s']),2)}
    #save cross-train scores for iteration
    param_acc_2[name] = avg_acc
    param_prec_2[name] = avg_prec
    param_rec_2[name] = avg_rec  
    param_f1_2[name] = avg_f1

    #%% Find best parameters for logistic regression
#scale the data using standard scaler
std_scaler = preprocessing.StandardScaler()
train_scaled = std_scaler.fit_transform(train.drop(['call'],axis=1))
### Split into train/test sets, and then cross-validation sets
Y = train.call.replace({'a':0,'f':1,'s':2}).astype(int)
#X = train.drop(['call'],axis=1).astype(float).values #get features
X = train_scaled
#fix random seed for reproducibility
seed = 999
np.random.seed(seed)
#split into train/test
x_train,x_test,y_train,y_test = train_test_split(X,Y,test_size=0.1,random_state=seed)
#set cross-validation params
kfold = KFold(n_splits=5,shuffle=True,random_state=seed)


''' default settings
penalty = l2
C = 1
solver = 'lbfgs' for l2 and 'saga' for l1/enet
'''

params = {'penalty':'l2','dual':False,'tol':0.0001,'C':1.0,'fit_intercept':True,'intercept_scaling':1,'class_weight':None,'random_state':None,'solver':'lbfgs',
          'max_iter':5000,'multi_class':'auto','verbose':0,'warm_start':False,'n_jobs':4,'l1_ratio':None}

#test parameters in k-fold cross-validation
param_acc_lr = {}
param_prec_lr = {}
param_rec_lr = {}
param_f1_lr = {}
for name,setting in [('lasso-0.001',{'penalty':'l1','C':0.001,'solver':'saga'}),
                     ('lasso-0.01',{'penalty':'l1','C':0.01,'solver':'saga'}),
                     ('lasso-0.1',{'penalty':'l1','C':0.1,'solver':'saga'}),
                     ('lasso-1',{'penalty':'l1','C':1,'solver':'saga'}),
                     ('lasso-10',{'penalty':'l1','C':10,'solver':'saga'}),
                     ('lasso-100',{'penalty':'l1','C':100,'solver':'saga'}),
                     ('ridge-0.001',{'penalty':'l2','C':0.001,'solver':'lbfgs'}),
                     ('ridge-0.01',{'penalty':'l2','C':0.01,'solver':'lbfgs'}),
                     ('ridge-0.1',{'penalty':'l2','C':0.1,'solver':'lbfgs'}),
                     ('ridge-1',{'penalty':'l2','C':1,'solver':'lbfgs'}),
                     ('ridge-10',{'penalty':'l2','C':10,'solver':'lbfgs'}),
                     ('ridge-100',{'penalty':'l2','C':100,'solver':'lbfgs'}),
                     ('elnet-0.001',{'penalty':'elasticnet','C':0.001,'solver':'saga'}),
                     ('elnet-0.01',{'penalty':'elasticnet','C':0.01,'solver':'saga'}),
                     ('elnet-0.1',{'penalty':'elasticnet','C':0.1,'solver':'saga'}),
                     ('elnet-1',{'penalty':'elasticnet','C':1,'solver':'saga'}),
                     ('elnet-10',{'penalty':'elasticnet','C':10,'solver':'saga'}),
                     ('elnet-100',{'penalty':'elasticnet','C':100,'solver':'saga'})]:

    #update the model parameters
    params.update(setting)
    print(name,params)

    #init model
    lrm = LogisticRegression(**params)

    #cross-validation
    cv_acc = []
    cv_prec = {'a':[],'f':[],'s':[]}
    cv_rec = {'a':[],'f':[],'s':[]}
    cv_f1 = {'a':[],'f':[],'s':[]}
    for trainix,testix in kfold.split(X=X,y=Y):
        
        #split into train,test sets
        x_train = X[trainix]
        y_train = Y.iloc[trainix]
        x_test = X[testix]
        y_test = Y.iloc[testix]
        
        #fit model and predict
        lrm.fit(x_train,y_train)
        pred_labels = lrm.predict(x_test)
        
        #assess performance
        acc = round(metrics.accuracy_score(y_test, pred_labels),2)
        ps = metrics.precision_score(y_test,pred_labels,average=None)
        rs = metrics.recall_score(y_test,pred_labels,average=None)
        fs = metrics.f1_score(y_test,pred_labels,average=None)
    
        #save values
        cv_acc.append(acc)
        cv_prec['a'] = cv_prec['a'] + [round(ps[0],2)]
        cv_prec['f'] = cv_prec['f'] + [round(ps[1],2)]
        cv_prec['s'] = cv_prec['s'] + [round(ps[2],2)]
        cv_rec['a'] = cv_rec['a'] + [round(rs[0],2)]
        cv_rec['f'] = cv_rec['f'] + [round(rs[1],2)]
        cv_rec['s'] = cv_rec['s'] + [round(rs[2],2)]
        cv_f1['a'] = cv_f1['a'] + [round(fs[0],2)]
        cv_f1['f'] = cv_f1['f'] + [round(fs[1],2)]
        cv_f1['s'] = cv_f1['s'] + [round(fs[2],2)]
        
    #get average cross-train scores
    avg_acc = round(np.mean(cv_acc),2)
    avg_prec = {'a':round(np.mean(cv_prec['a']),2),'f':round(np.mean(cv_prec['f']),2),'s':round(np.mean(cv_prec['s']),2)}
    avg_rec = {'a':round(np.mean(cv_rec['a']),2),'f':round(np.mean(cv_rec['f']),2),'s':round(np.mean(cv_rec['s']),2)}
    avg_f1 = {'a':round(np.mean(cv_f1['a']),2),'f':round(np.mean(cv_f1['f']),2),'s':round(np.mean(cv_f1['s']),2)}
    #save cross-train scores for iteration
    param_acc_lr[name] = avg_acc
    param_prec_lr[name] = avg_prec
    param_rec_lr[name] = avg_rec  
    param_f1_lr[name] = avg_f1
    
#test parameters in elnet k-fold cross-validation
params = {'penalty':'l2','dual':False,'tol':0.0001,'C':1.0,'fit_intercept':True,'intercept_scaling':1,'class_weight':None,'random_state':None,'solver':'lbfgs',
      'max_iter':5000,'multi_class':'auto','verbose':0,'warm_start':False,'n_jobs':4,'l1_ratio':None}
param_acc_lr_enet = {}
param_prec_lr_enet = {}
param_rec_lr_enet = {}
param_f1_lr_enet = {}
for name,setting in [('elnet-0.01-0.25',{'penalty':'elasticnet','C':0.01,'solver':'saga','l1_ratio':0.25}),
                     ('elnet-0.01-0.50',{'penalty':'elasticnet','C':0.01,'solver':'saga','l1_ratio':0.5}),
                     ('elnet-0.01-0.75',{'penalty':'elasticnet','C':0.01,'solver':'saga','l1_ratio':0.75}),
                     ('elnet-0.1-0.25',{'penalty':'elasticnet','C':0.1,'solver':'saga','l1_ratio':0.25}),
                     ('elnet-0.1-0.5',{'penalty':'elasticnet','C':0.1,'solver':'saga','l1_ratio':0.5}),
                     ('elnet-0.1-0.75',{'penalty':'elasticnet','C':0.1,'solver':'saga','l1_ratio':0.75}),
                     ('elnet-1-0.25',{'penalty':'elasticnet','C':1,'solver':'saga','l1_ratio':0.25}),
                     ('elnet-1-0.5',{'penalty':'elasticnet','C':1,'solver':'saga','l1_ratio':0.5}),
                     ('elnet-1-0.75',{'penalty':'elasticnet','C':1,'solver':'saga','l1_ratio':0.75}),
                     ('elnet-10-0.25',{'penalty':'elasticnet','C':10,'solver':'saga','l1_ratio':0.25}),
                     ('elnet-10-0.5',{'penalty':'elasticnet','C':10,'solver':'saga','l1_ratio':0.5}),
                     ('elnet-10-0.75',{'penalty':'elasticnet','C':10,'solver':'saga','l1_ratio':0.75}),
                    ]:

    params.update(setting)
    print(name,params)

    #init model
    lrm = LogisticRegression(**params)

    #cross-validation
    cv_acc = []
    cv_prec = {'a':[],'f':[],'s':[]}
    cv_rec = {'a':[],'f':[],'s':[]}
    cv_f1 = {'a':[],'f':[],'s':[]}
    for trainix,testix in kfold.split(X=X,y=Y):
        
        #split into train,test sets
        x_train = X[trainix]
        y_train = Y.iloc[trainix]
        x_test = X[testix]
        y_test = Y.iloc[testix]
        
        #fit model and predict
        lrm.fit(x_train,y_train)
        pred_labels = lrm.predict(x_test)
        
        #assess performance
        acc = round(metrics.accuracy_score(y_test, pred_labels),2)
        ps = metrics.precision_score(y_test,pred_labels,average=None)
        rs = metrics.recall_score(y_test,pred_labels,average=None)
        fs = metrics.f1_score(y_test,pred_labels,average=None)
    
        #save values
        cv_acc.append(acc)
        cv_prec['a'] = cv_prec['a'] + [round(ps[0],2)]
        cv_prec['f'] = cv_prec['f'] + [round(ps[1],2)]
        cv_prec['s'] = cv_prec['s'] + [round(ps[2],2)]
        cv_rec['a'] = cv_rec['a'] + [round(rs[0],2)]
        cv_rec['f'] = cv_rec['f'] + [round(rs[1],2)]
        cv_rec['s'] = cv_rec['s'] + [round(rs[2],2)]
        cv_f1['a'] = cv_f1['a'] + [round(fs[0],2)]
        cv_f1['f'] = cv_f1['f'] + [round(fs[1],2)]
        cv_f1['s'] = cv_f1['s'] + [round(fs[2],2)]
        
    #get average cross-train scores
    avg_acc = round(np.mean(cv_acc),2)
    avg_prec = {'a':round(np.mean(cv_prec['a']),2),'f':round(np.mean(cv_prec['f']),2),'s':round(np.mean(cv_prec['s']),2)}
    avg_rec = {'a':round(np.mean(cv_rec['a']),2),'f':round(np.mean(cv_rec['f']),2),'s':round(np.mean(cv_rec['s']),2)}
    avg_f1 = {'a':round(np.mean(cv_f1['a']),2),'f':round(np.mean(cv_f1['f']),2),'s':round(np.mean(cv_f1['s']),2)}
    #save cross-train scores for iteration
    param_acc_lr_enet[name] = avg_acc
    param_prec_lr_enet[name] = avg_prec
    param_rec_lr_enet[name] = avg_rec  
    param_f1_lr_enet[name] = avg_f1

    ### results of model cross-validations

#RF
'''
Use 500 trees
cv_acc: 0.9
cv_prec: {'a': 0.86, 'f': 0.92, 's': 0.9}
cv_recall: {'a': 0.79, 'f': 0.92, 's': 0.94}
cv_f1: {'a': 0.83, 'f': 0.92, 's': 0.92}
'''

#RBF-SVM
'''
use kernel parameter of 10
--- cv avg accuracy
- 0.1: 0.84,
- 1: 0.87,
- 10: 0.87,
- 100: 0.87}
--- cv avg precision
- 0.1: {'a': 0.78, 'f': 0.89, 's': 0.85},
- 1: {'a': 0.8, 'f': 0.92, 's': 0.87},
- 10: {'a': 0.82, 'f': 0.91, 's': 0.88},
- 100: {'a': 0.8, 'f': 0.91, 's': 0.88}}
--- cv avg recall
- 0.1: {'a': 0.7, 'f': 0.9, 's': 0.9},
- 1: {'a': 0.74, 'f': 0.92, 's': 0.91},
- 10: {'a': 0.74, 'f': 0.92, 's': 0.93},
- 100: {'a': 0.75, 'f': 0.9, 's': 0.92}}
--- cv avg f1 score
- 0.1: {'a': 0.74, 'f': 0.9, 's': 0.87},
- 1: {'a': 0.77, 'f': 0.92, 's': 0.89},
- 10: {'a': 0.78, 'f': 0.91, 's': 0.9},
- 100: {'a': 0.78, 'f': 0.9, 's': 0.9}}
'''

#gradient boosted trees
'''
best params:
'0.1,1000,6': 0.9} --> 0.9 accuracy
---accuracy
{'0.1,250,4': 0.89,
 '0.1,500,4': 0.89,
 '0.1,1000,4': 0.89,
 '0.1,250,5': 0.89,
 '0.1,500,5': 0.89,
 '0.1,1000,5': 0.89,
 '0.1,250,6': 0.89,
 '0.1,500,6': 0.89,
 '0.1,1000,6': 0.9}
---precision
{'0.1,250,4': {'a': 0.84, 'f': 0.92, 's': 0.9},
 '0.1,500,4': {'a': 0.84, 'f': 0.92, 's': 0.9},
 '0.1,1000,4': {'a': 0.84, 'f': 0.92, 's': 0.9},
 '0.1,250,5': {'a': 0.85, 'f': 0.92, 's': 0.9},
 '0.1,500,5': {'a': 0.85, 'f': 0.92, 's': 0.9},
 '0.1,1000,5': {'a': 0.84, 'f': 0.92, 's': 0.91},
 '0.1,250,6': {'a': 0.85, 'f': 0.92, 's': 0.91},
 '0.1,500,6': {'a': 0.85, 'f': 0.92, 's': 0.91},
 '0.1,1000,6': {'a': 0.85, 'f': 0.93, 's': 0.91}}
---recall
{'0.1,250,4': {'a': 0.79, 'f': 0.92, 's': 0.93},
 '0.1,500,4': {'a': 0.8, 'f': 0.92, 's': 0.93},
 '0.1,1000,4': {'a': 0.8, 'f': 0.92, 's': 0.93},
 '0.1,250,5': {'a': 0.79, 'f': 0.92, 's': 0.94},
 '0.1,500,5': {'a': 0.8, 'f': 0.92, 's': 0.93},
 '0.1,1000,5': {'a': 0.8, 'f': 0.92, 's': 0.93},
 '0.1,250,6': {'a': 0.8, 'f': 0.92, 's': 0.94},
 '0.1,500,6': {'a': 0.8, 'f': 0.93, 's': 0.93},
 '0.1,1000,6': {'a': 0.8, 'f': 0.92, 's': 0.94}}
---f1
{'0.1,250,4': {'a': 0.82, 'f': 0.92, 's': 0.92},
 '0.1,500,4': {'a': 0.82, 'f': 0.92, 's': 0.92},
 '0.1,1000,4': {'a': 0.82, 'f': 0.92, 's': 0.92},
 '0.1,250,5': {'a': 0.82, 'f': 0.92, 's': 0.92},
 '0.1,500,5': {'a': 0.82, 'f': 0.92, 's': 0.92},
 '0.1,1000,5': {'a': 0.82, 'f': 0.92, 's': 0.92},
 '0.1,250,6': {'a': 0.82, 'f': 0.92, 's': 0.92},
 '0.1,500,6': {'a': 0.82, 'f': 0.93, 's': 0.92},
 '0.1,1000,6': {'a': 0.82, 'f': 0.93, 's': 0.92}}
'''

#Logistic regression
'''
Logistic regression results
------------Accuracy------------
 'elnet-0.01-0.25': 0.83,
 'elnet-0.01-0.50': 0.83,
 'elnet-0.01-0.75': 0.82,
 'elnet-0.1-0.25': 0.84,
 'elnet-0.1-0.5': 0.84,
 'elnet-0.1-0.75': 0.84,
 'elnet-1-0.25': 0.84,
 'elnet-1-0.5': 0.84,
 'elnet-1-0.75': 0.84,
 'elnet-10-0.25': 0.84,
 'elnet-10-0.5': 0.84,
 'elnet-10-0.75': 0.84
 'lasso-0.001': 0.81,
 'lasso-0.01': 0.82,
 'lasso-0.1': 0.84,
 'lasso-1': 0.84,
 'lasso-10': 0.84,
 'lasso-100': 0.84,
 'ridge-0.001': 0.82,
 'ridge-0.01': 0.83,
 'ridge-0.1': 0.84,
 'ridge-1': 0.84,
 'ridge-10': 0.84,
 'ridge-100': 0.84
------------Precision------------
{'elnet-0.01-0.25': {'a': 0.75, 'f': 0.91, 's': 0.83},
 'elnet-0.01-0.50': {'a': 0.74, 'f': 0.91, 's': 0.82},
 'elnet-0.01-0.75': {'a': 0.74, 'f': 0.91, 's': 0.82},
 'elnet-0.1-0.25': {'a': 0.76, 'f': 0.91, 's': 0.84},
 'elnet-0.1-0.5': {'a': 0.76, 'f': 0.91, 's': 0.84},
 'elnet-0.1-0.75': {'a': 0.76, 'f': 0.91, 's': 0.84},
 'elnet-1-0.25': {'a': 0.76, 'f': 0.91, 's': 0.85},
 'elnet-1-0.5': {'a': 0.76, 'f': 0.91, 's': 0.85},
 'elnet-1-0.75': {'a': 0.76, 'f': 0.91, 's': 0.85},
 'elnet-10-0.25': {'a': 0.76, 'f': 0.91, 's': 0.85},
 'elnet-10-0.5': {'a': 0.76, 'f': 0.91, 's': 0.85},
 'elnet-10-0.75': {'a': 0.76, 'f': 0.91, 's': 0.85}}
{'lasso-0.001': {'a': 0.81, 'f': 0.89, 's': 0.77},
 'lasso-0.01': {'a': 0.74, 'f': 0.91, 's': 0.82},
 'lasso-0.1': {'a': 0.76, 'f': 0.91, 's': 0.84},
 'lasso-1': {'a': 0.76, 'f': 0.91, 's': 0.85},
 'lasso-10': {'a': 0.76, 'f': 0.91, 's': 0.85},
 'lasso-100': {'a': 0.76, 'f': 0.91, 's': 0.85},
 'ridge-0.001': {'a': 0.77, 'f': 0.91, 's': 0.81},
 'ridge-0.01': {'a': 0.75, 'f': 0.91, 's': 0.83},
 'ridge-0.1': {'a': 0.75, 'f': 0.91, 's': 0.84},
 'ridge-1': {'a': 0.76, 'f': 0.91, 's': 0.85},
 'ridge-10': {'a': 0.76, 'f': 0.91, 's': 0.85},
 'ridge-100': {'a': 0.76, 'f': 0.91, 's': 0.85}}
------------Recall------------
{'elnet-0.01-0.25': {'a': 0.66, 'f': 0.88, 's': 0.9},
 'elnet-0.01-0.50': {'a': 0.65, 'f': 0.88, 's': 0.9},
 'elnet-0.01-0.75': {'a': 0.65, 'f': 0.87, 's': 0.89},
 'elnet-0.1-0.25': {'a': 0.68, 'f': 0.88, 's': 0.9},
 'elnet-0.1-0.5': {'a': 0.68, 'f': 0.88, 's': 0.9},
 'elnet-0.1-0.75': {'a': 0.68, 'f': 0.88, 's': 0.9},
 'elnet-1-0.25': {'a': 0.69, 'f': 0.88, 's': 0.9},
 'elnet-1-0.5': {'a': 0.69, 'f': 0.88, 's': 0.9},
 'elnet-1-0.75': {'a': 0.69, 'f': 0.88, 's': 0.9},
 'elnet-10-0.25': {'a': 0.69, 'f': 0.88, 's': 0.9},
 'elnet-10-0.5': {'a': 0.69, 'f': 0.88, 's': 0.9},
 'elnet-10-0.75': {'a': 0.69, 'f': 0.88, 's': 0.9}}
{'lasso-0.001': {'a': 0.54, 'f': 0.83, 's': 0.94},
 'lasso-0.01': {'a': 0.65, 'f': 0.87, 's': 0.89},
 'lasso-0.1': {'a': 0.68, 'f': 0.88, 's': 0.9},
 'lasso-1': {'a': 0.69, 'f': 0.88, 's': 0.9},
 'lasso-10': {'a': 0.69, 'f': 0.88, 's': 0.9},
 'lasso-100': {'a': 0.69, 'f': 0.88, 's': 0.9},
 'ridge-0.001': {'a': 0.61, 'f': 0.87, 's': 0.91},
 'ridge-0.01': {'a': 0.66, 'f': 0.88, 's': 0.9},
 'ridge-0.1': {'a': 0.68, 'f': 0.88, 's': 0.9},
 'ridge-1': {'a': 0.69, 'f': 0.88, 's': 0.9},
 'ridge-10': {'a': 0.69, 'f': 0.88, 's': 0.91},
 'ridge-100': {'a': 0.69, 'f': 0.88, 's': 0.91}}
------------F1 score------------
{'elnet-0.01-0.25': {'a': 0.7, 'f': 0.9, 's': 0.87},
 'elnet-0.01-0.50': {'a': 0.69, 'f': 0.9, 's': 0.86},
 'elnet-0.01-0.75': {'a': 0.69, 'f': 0.89, 's': 0.86},
 'elnet-0.1-0.25': {'a': 0.72, 'f': 0.9, 's': 0.87},
 'elnet-0.1-0.5': {'a': 0.72, 'f': 0.9, 's': 0.87},
 'elnet-0.1-0.75': {'a': 0.72, 'f': 0.9, 's': 0.87},
 'elnet-1-0.25': {'a': 0.72, 'f': 0.9, 's': 0.88},
 'elnet-1-0.5': {'a': 0.72, 'f': 0.9, 's': 0.88},
 'elnet-1-0.75': {'a': 0.72, 'f': 0.9, 's': 0.88},
 'elnet-10-0.25': {'a': 0.72, 'f': 0.9, 's': 0.88},
 'elnet-10-0.5': {'a': 0.72, 'f': 0.9, 's': 0.88},
 'elnet-10-0.75': {'a': 0.72, 'f': 0.9, 's': 0.88}}
{'lasso-0.001': {'a': 0.65, 'f': 0.86, 's': 0.85},
 'lasso-0.01': {'a': 0.69, 'f': 0.89, 's': 0.86},
 'lasso-0.1': {'a': 0.72, 'f': 0.9, 's': 0.87},
 'lasso-1': {'a': 0.72, 'f': 0.9, 's': 0.88},
 'lasso-10': {'a': 0.72, 'f': 0.9, 's': 0.88},
 'lasso-100': {'a': 0.72, 'f': 0.9, 's': 0.88},
 'ridge-0.001': {'a': 0.68, 'f': 0.89, 's': 0.86},
 'ridge-0.01': {'a': 0.7, 'f': 0.9, 's': 0.87},
 'ridge-0.1': {'a': 0.72, 'f': 0.9, 's': 0.87},
 'ridge-1': {'a': 0.72, 'f': 0.9, 's': 0.88},
 'ridge-10': {'a': 0.72, 'f': 0.9, 's': 0.88},
 'ridge-100': {'a': 0.73, 'f': 0.9, 's': 0.88}}
'''
