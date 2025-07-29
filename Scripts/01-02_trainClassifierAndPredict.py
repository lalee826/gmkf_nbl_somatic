### Project description
'''
- This script trains a new neural network on DeepSVR data using only features that can be obtained for GMKF data
- We remove liquid tumor variants and collapse the feature set to 50 extractable features from our own data set
- The deep learning package is Keras
'''

### Import packages
import pandas as pd
import numpy as np
import os
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import train_test_split
from sklearn import preprocessing as preprocessing
from sklearn import metrics as metrics
from tensorflow import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
from tensorflow.keras.regularizers import l2
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier

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
#set cross-val params
seed = 999
kfold = KFold(n_splits=5,shuffle=True,random_state=np.random.seed(seed))

#test parameters of nn structure
nnodes = [25,50]
nlayers = [2,3]
activations = ['relu','sigmoid','tanh']

#store cross-validation metrics
nnCV = {}

#change y train array from one-hot encoded to single dimension
y_train_1d = np.argmax(y_train, axis=1)

#get feature number
nFeatures = x_train.shape[1]

for act in activations:
    for layer in nlayers:
        for nnum in nnodes:
            print('activation function: {}'.format(act))
            print('layer number: {}'.format(layer))
            print('nodes per layer: {}'.format(nnum))
            
            #name the parameter set
            name = 'activation: ' + act + ', num_layers: ' + str(layer) + ', num_nodes: ' + str(nnum)
            
            #build the model using training params
            def buildModelCV(actf,nlayer,nodes):
                model = keras.models.Sequential()
                model.add(Dense(50,input_dim=nFeatures,kernel_initializer=keras.initializers.RandomUniform(minval=-0.05,maxval=0.05,seed=seed),
                            activation=actf,kernel_regularizer=keras.regularizers.l2(0.001))) #input and output of initial layer are both matrices of dimensions [*,50]
                for i in range(nlayer): #add layers
                    model.add(Dense(nodes,activation=actf,kernel_regularizer=l2(0.001)))
                model.add(Dense(3,kernel_initializer='random_uniform',activation='softmax')) #add a final softmax layer corresponding to three labels
                model.compile(loss='sparse_categorical_crossentropy', optimizer='adam', metrics=['accuracy']) #compile the model
                return(model)
            
            #store cv scores
            acc = []
            prec = {0:[],1:[],2:[]}
            recall = {0:[],1:[],2:[]}
            f1 = {0:[],1:[],2:[]}
            
            #test parameters in cross-validation
            for trainix,testix in kfold.split(x_train):
                
                #split cross val data
                xtrain,ytrain = x_train[trainix,:], y_train_1d[trainix]
                xtest,ytest = x_train[testix,:], y_train_1d[testix]           
                
                #build nn with specific params
                nnModel = KerasClassifier(build_fn=lambda: buildModelCV(act,layer,nnum), epochs=100, batch_size=1500, verbose=0)
                
                #fit the model of cross val train data
                fitmodel = nnModel.fit(xtrain, ytrain)
                      
                #evaluate model on test set
                predicted = nnModel.predict(xtest, verbose=0)
                      
                #classify test set
                #predicted = np.array([list(a).index(max(list(a))) for a in list(pred_labels)])
                
                #score the predictions
                label_binarizer = preprocessing.LabelBinarizer()
                label_binarizer.fit(range(max(predicted)+1))
                ypredOneHot = label_binarizer.transform(predicted)
                ytestOneHot = label_binarizer.transform(ytest)
                
                #evaluate performance
                ac = round(metrics.accuracy_score(ytestOneHot,ypredOneHot),2)
                pr = metrics.precision_score(ytestOneHot,ypredOneHot,average=None).round(2)
                re = metrics.recall_score(ytestOneHot,ypredOneHot,average=None).round(2)
                f1s = metrics.f1_score(ytestOneHot,ypredOneHot,average=None).round(2)
        
                
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
            nnCV[name] = {'name':name,
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
                
### train our neural network with best parameters
'''
best parameters in cross-validation:
 'activation: relu, num_layers: 3, num_nodes: 50': {'avgAccuracy': 0.85,
  'avgF1_a': 0.9,
  'avgF1_f': 0.73,
  'avgF1_s': 0.88,
  'avgPrecision_a': 0.91,
  'avgPrecision_f': 0.78,
  'avgPrecision_s': 0.86,
  'avgRecall_a': 0.9,
  'avgRecall_f': 0.69,
  'avgRecall_s': 0.91,
'''
#configure the model based on the data
n_input, n_classes = x_train.shape[1], y_train.shape[1]
#define model
model = keras.models.Sequential() #connected feed forward
#add first layer 50 nodes wide
model.add(Dense(50, input_dim=n_input, activation='relu', kernel_initializer='normal',kernel_regularizer=keras.regularizers.l2(0.001)))
#add three layers 25 nodes wide each
for i in range(3): #add layers
    model.add(Dense(50,activation='relu',kernel_regularizer=l2(0.001)))
#add 3 layers 6 nodes wide each
#for i in range(10): #add layers
#    model.add(Dense(6,activation='sigmoid',kernel_regularizer=l2(0.001)))
#add a final softmax layer
model.add(Dense(n_classes,kernel_initializer='normal',activation='softmax'))
#determine how to compile model
model.compile(loss='sparse_categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
# fit model on train set
model.fit(x_train,y_train_1d, epochs=100, verbose=0)#, batch_size=2000)

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
#min max scale the data
mm_scaler = preprocessing.MinMaxScaler()
data_scaled = mm_scaler.fit_transform(nblData)
data_scaled = pd.DataFrame(data_scaled, index=nblData.index,
                      columns=nblData.columns)

#make predictions using trained model
predLabelsNBL = model.predict(data_scaled, verbose=0)

#use max softmax probability to assign labels
predicted = np.array([list(a).index(max(list(a))) for a in list(predLabelsNBL)])
#look at prediction counts
predCounts = np.bincount(predicted)

#add classifier predictions to each variant
nblData['prediction'] = predicted
nblData['prediction'].replace({0:'ambiguous',1:'fail',2:'somatic'},inplace=True)

#join back to consensus and unique variant mutations files
consensusVars = consensusVars.merge(nblData,on='mutid',how='left')
uniqueVars = uniqueVars.merge(nblData,on='mutid',how='left')
