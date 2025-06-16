#%% Project description
'''
- This script trains a new neural network on DeepSVR data using only features that can be obtained for GMKF data
- We remove liquid tumor variants and collapse the feature set to 50 extractable features from our own data set
- The deep learning package is Keras

Author: Alex Lee
Date: 9/19/2019
'''

#%%Import packages
import pandas as pd
import numpy as np
import os
import sklearn
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
from sklearn import metrics
import matplotlib.pyplot as plt
from itertools import cycle
from tensorflow import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier
from tensorflow.keras.regularizers import l2

workdir = '/rocker-build/gmkf_nbl_somatic/'

#%% Import and prepare training data
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

#import our own data
nblData = workdir + 'Data/gmkf_matched_features.pkl'
tmbdata = pd.read_pickle(nblData)
tmbdata = tmbdata.rename(columns={'normal_ref_avg_num_mismaches_as_fraction':'normal_ref_avg_num_mismatches_as_fraction',
                      'tumor_ref_avg_num_mismaches_as_fraction':'tumor_ref_avg_num_mismatches_as_fraction',
                      'normal_var_avg_num_mismaches_as_fraction':'normal_var_avg_num_mismatches_as_fraction',
                      'tumor_var_avg_num_mismaches_as_fraction':'tumor_var_avg_num_mismatches_as_fraction'})
tmbdata = tmbdata[trainCols[1:].tolist()]

mm_scaler = preprocessing.MinMaxScaler()
data_scaled = mm_scaler.fit_transform(tmbdata)
data_scaled = pd.DataFrame(data_scaled, index=tmbdata.index,
                      columns=tmbdata.columns)

trainss = train.iloc[:100,:]
datass = data_scaled.iloc[:100,:]

train.shape #(32308, 59)
data_scaled.shape #(12481, 58)

#%% Split train data into features and labels
train_y = pd.get_dummies(train.call).astype(float).values #get labels
'''
0 - ambiguous
1 - fail
2 - somatic
'''
train_x = train.drop(['call'],axis=1).astype(float).values #get features
#fix random seed for reproducibility
seed = 999
np.random.seed(seed)

#%% train our neural network
#we only want simple feed-forward layers (type Dense)
#configure the model based on the data
n_input, n_classes = train_x.shape[1], train_y.shape[1]
#define model
model = keras.models.Sequential() #connected feed forward
#add first layer 50 nodes wide
model.add(Dense(50, input_dim=n_input, activation='relu', kernel_initializer='normal',kernel_regularizer=keras.regularizers.l2(0.001)))
#add three layers 25 nodes wide each
for i in range(3): #add layers
    model.add(Dense(25,activation='relu',kernel_regularizer=l2(0.001)))
#add 3 layers 6 nodes wide each
#for i in range(10): #add layers
#    model.add(Dense(6,activation='sigmoid',kernel_regularizer=l2(0.001)))
#add a final softmax layer
model.add(Dense(n_classes,kernel_initializer='normal',activation='softmax'))
#determine how to compile model
model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
# fit model on train set
fitmodel = model.fit(train_x,train_y, epochs=100, verbose=1)#, batch_size=2000)

#%% Make predictions
#evaluate model on test set
data_pred_labels = model.predict(data_scaled.iloc[:,:58], verbose=1)
predicted = np.array([list(a).index(max(list(a))) for a in list(data_pred_labels)])
#look at prediction counts
print(np.bincount(predicted))
'''
counts are shown as ambiguous, fail, and somatic
'''

#add classifier predictions to each variant
data_scaled['classification'] = predicted

#### CROSS VALIDATION OF BEST MODEL PARAMETERS ####
#%% Split into train/test sets, and then cross-validation sets
Y = pd.get_dummies(train.call).astype(float).values #get labels
X = train.drop(['call'],axis=1).astype(float).values #get features
#fix random seed for reproducibility
seed = 13
np.random.seed(seed)
#split into train/test
x_train,x_test,y_train,y_test = sklearn.model_selection.train_test_split(X,Y,test_size=0.15,random_state=seed)

#%% train neural network and evaluate performance
#usesimple feed-forward layers (type Dense)

#configure the model based on the data
n_input, n_classes = x_train.shape[1], y_test.shape[1]
#define model
model = keras.models.Sequential() #connected feed forward
#add first layer 50 nodes wide
model.add(Dense(50, input_dim=n_input, activation='relu', kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.001)))
#add three layers 25 nodes wide each
for i in range(3): #add layers
    model.add(Dense(25,activation='relu',kernel_regularizer=l2(0.001)))
#add 3 layers 6 nodes wide each
#for i in range(10): #add layers
#    model.add(Dense(6,activation='sigmoid',kernel_regularizer=l2(0.001)))
#add a final softmax layer
model.add(Dense(n_classes,kernel_initializer='random_uniform',activation='softmax'))
#determine how to compile model
model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
# fit model on train set
fitmodel = model.fit(x_train, y_train, epochs=100, verbose=1)
#evaluate model on test set
_, test_acc = model.evaluate(x_test, y_test, verbose=1)
pred_labels = model.predict(x_test, verbose=1)
#%%classify test set
predicted = np.array([list(a).index(max(list(a))) for a in list(pred_labels)])
#score the predictions
label_binarizer = preprocessing.LabelBinarizer()
label_binarizer.fit(range(max(predicted)+1))
predicted_transformed = label_binarizer.transform(predicted)

#%%check results
print('Model accuracy:\t' +  
      str(sklearn.metrics.accuracy_score(y_test, predicted_transformed)) + 
      '\nModel Performance report\n' + 
      str(sklearn.metrics.classification_report(y_test, predicted_transformed)))

#%% set cross-val params
kfold = sklearn.model_selection.KFold(n_splits=10,shuffle=True,random_state=seed)
nnodes, nlayers, activations = [25,50], [1,2,3], ['relu','sigmoid','tanh']
nnodes, nlayers, activations = [50], [3], ['relu']
#histories = {}
#test_results = {}
cv_scores = {}
probabilities_d = {}
#change y train array from one-hot encoded to single dimension
y_train_1d = np.argmax(y_train, axis=1)
input_dims = x_train.shape[1]

for act in activations:
    for layer in nlayers:
        for nnum in nnodes:
            print('activation function: {}'.format(act))
            print('layer number: {}'.format(layer))
            print('nodes per layer: {}'.format(nnum))
            #set training params
            def buildModelCV():
                model = keras.models.Sequential()
                model.add(Dense(50,input_dim=input_dims,kernel_initializer=keras.initializers.RandomUniform(minval=-0.05,maxval=0.05,seed=seed),
                            activation='relu',kernel_regularizer=keras.regularizers.l2(0.001))) #input and output of initial layer are both matrices of dimensions [*,50]
                for i in range(layer): #add layers
                    model.add(Dense(nnum,activation=act,kernel_regularizer=l2(0.001)))
                model.add(Dense(3,kernel_initializer='random_uniform',activation='softmax')) #add a final softmax layer corresponding to three labels
                model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy']) #compile the model
                return(model)
            nnEstimator = KerasClassifier(build_fn=buildModelCV, epochs=100, batch_size=1500, verbose=0)
            probabilities = cross_val_predict(nnEstimator, x_train, y_train_1d, cv=kfold, method='predict_proba')
            #score the predictions
            predicted = np.array([list(a).index(max(list(a))) for a in list(probabilities)])
            label_binarizer = preprocessing.LabelBinarizer()
            label_binarizer.fit(range(max(predicted)+1))
            predicted_transformed = label_binarizer.transform(predicted)
            cv_acc = sklearn.metrics.accuracy_score(y_train, predicted_transformed)
            probabilities_d.update({act + '_' + str(nnum) + '_' + str(layer) + 'layers':probabilities})
            cv_scores.update({act + '_' + str(nnum) + '_' + str(layer) + 'layers':cv_acc})
            print('number of layers = %d; number of nodes = %d: %.3f' % (layer, nnum, cv_acc))

#%%top performance is in relu/50/3
probabilities = probabilities_d['relu_50_3layers']
predicted = np.array([list(a).index(max(list(a))) for a in list(probabilities)])
#score the predictions
label_binarizer = preprocessing.LabelBinarizer()
label_binarizer.fit(range(max(predicted)+1))
predicted_transformed = label_binarizer.transform(predicted)

#%% evaluate performance of best model
print('Cross validation accuracy:')
print('\t', sklearn.metrics.accuracy_score(y_train, predicted_transformed))
print('\nCross validation classification report\n')
print(sklearn.metrics.classification_report(y_train, predicted_transformed))

#ROC for three classes using RELU/3/50
def create_roc_curve(Y, probabilities, class_lookup, title, ax):
    '''Create ROC curve to compare multiclass model performance.

    Parameters:
        Y (numpy.array): Truth labels
        probabilities (numpy.array): Output of model for each class
        class_lookup (dict): lookup hash of truth labels
        title (str): Plot title
    '''
    n_classes = Y.shape[1]
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    ax.set_title(title)
    if n_classes == 3:
        colors = cycle(['orange', 'red', 'black'])
    else:
        colors = cycle(['orange', 'red', 'aqua', 'black'])
    for i, color in zip(range(n_classes), colors):
        fpr[i], tpr[i], _ = metrics.roc_curve(Y[:, i], probabilities[:, i])
        roc_auc[i] = metrics.auc(fpr[i], tpr[i])
        ax.plot(fpr[i], tpr[i], color=color,
                label='ROC curve of class {0} (area = {1:0.2f})'.format(
                    class_lookup[i], roc_auc[i]))
    ax.plot([0, 1], [0, 1], 'k--')
    ax.set_xlim([0, 1.0])
    ax.set_ylim([0.0, 1.0])
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.legend(loc="lower right")
    plt.show()

class_lookup = {0: 'Ambiguous', 1: 'Fail', 2: 'Somatic'}
fig, ax = plt.subplots()
create_roc_curve(y_train, probabilities, class_lookup, 'Three Class Reciever '
                'Operating Characteristic Curve',ax)


#output the final count of labels on our test data set
predCounts = data_scaled.groupby('classification').size()
