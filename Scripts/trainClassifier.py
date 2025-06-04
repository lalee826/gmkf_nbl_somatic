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






































                            ##############################################
                            ''' TRAINING AND CROSS-VALIDATION BELOW '''
                            ##############################################




#%% Split into train/test sets, and then cross-validation sets
Y = pd.get_dummies(train.call).astype(float).values #get labels
X = train.drop(['call'],axis=1).astype(float).values #get features
#fix random seed for reproducibility
seed = 13
np.random.seed(seed)
#split into train/test
x_train,x_test,y_train,y_test = sklearn.model_selection.train_test_split(X,Y,test_size=0.15,random_state=seed)
#set cross-validation params
#kfold = sklearn.model_selection.KFold(n_splits=10,shuffle=True,random_state=seed)

#%% train our neural network and evaluate performance
#we only want simple feed-forward layers (type Dense)

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
fitmodel = model.fit(x_train, y_train, epochs=1000, verbose=1)
#evaluate model on test set
_, test_acc = model.evaluate(x_test, y_test, verbose=1)
pred_labels = model.predict(x_test, verbose=1)
#%%classify test set
predicted = np.array([list(a).index(max(list(a))) for a in list(pred_labels)])
#score the predictions
label_binarizer = preprocessing.LabelBinarizer()
label_binarizer.fit(range(max(predicted)+1))
predicted_transformed = label_binarizer.transform(predicted)


#%%
print('Model accuracy:')
print('\t', sklearn.metrics.accuracy_score(y_test, predicted_transformed))
print('\nModel Performance report\n')
print(sklearn.metrics.classification_report(y_test, predicted_transformed))

#save the trained model
# save model and architecture to single file
model.save("/Volumes/SSD/Alex/model.tmb")


'''
def buildModelTrain():
    model = keras.models.Sequential()
    model.add(Dense(50,input_dim=input_dims,kernel_initializer=keras.initializers.RandomUniform(minval=-0.05,maxval=0.05,seed=seed),
                activation='relu',kernel_regularizer=keras.regularizers.l2(0.001))) #input and output of initial layer are both matrices of dimensions [*,50]
    for i in range(3): #add layers
        model.add(Dense(25,activation='relu',kernel_regularizer=l2(0.001)))
    model.add(Dense(3,kernel_initializer='random_uniform',activation='softmax')) #add a final softmax layer corresponding to three labels
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy']) #compile the model
    return(model)

nnEstimator = KerasClassifier(build_fn=buildModelTrain, epochs=1000, batch_size=1500, verbose=1)

def buildModel():
    model = keras.models.Sequential()
    model.add(Dense(50,input_dim=50,kernel_initializer=keras.initializers.RandomUniform(minval=-0.05,maxval=0.05,seed=seed),
                activation='relu',kernel_regularizer=keras.regularizers.l2(0.001))) #input and output of initial layer are both matrices of dimensions [*,50]
    for i in range(3): #add three layers, each with output of dimension [*,20] and relu activation fxn
        model.add(Dense(20,activation='relu',kernel_regularizer=l2(0.001)))
    model.add(Dense(3,kernel_initializer='random_uniform',activation='softmax')) #add a final softmax layer corresponding to three labels
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy']) #compile the model
    return(model)


def buildModelCV(nnodes,actf,n_input):
    print(nnodes,actf,n_input)
    model = keras.models.Sequential()
    model.add(Dense(50,input_dim=50,kernel_initializer=keras.initializers.RandomUniform(minval=-0.05,maxval=0.05,seed=seed),
                activation='relu',kernel_regularizer=keras.regularizers.l2(0.001))) #input and output of initial layer are both matrices of dimensions [*,50]
    model.add(Dense(3,kernel_initializer='random_uniform',activation='softmax')) #add a final softmax layer corresponding to three labels
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy']) #compile the model
    return(model)
'''


#set cross-val params
nnodes, nlayers, activations = [25,50], [1,2,3], ['relu','sigmoid','tanh']
nnodes, nlayers, activations = [50], [3], ['relu']
#histories = {}
#test_results = {}
cv_scores = {}
probabilities_d = {}

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
            nnEstimator = KerasClassifier(build_fn=buildModelCV, epochs=1000, batch_size=1500, verbose=0)
            #nnEstimator = KerasClassifier(build_fn=buildModel, epochs=1, batch_size=1500, verbose=1)
            #fitmodel, result, cvSC = evaluate_first_layer(nnum,act,x_train,y_train,x_test,y_test)
            #histories.update({act + '_' + str(nnum):fitmodel.history})
            #test_results.update({act + '_' + str(nnum):result})
            probabilities = cross_val_predict(nnEstimator, x_train, y_train, cv=kfold, method='predict_proba')
            #score the predictions
            predicted = np.array([list(a).index(max(list(a))) for a in list(probabilities)])
            label_binarizer = preprocessing.LabelBinarizer()
            label_binarizer.fit(range(max(predicted)+1))
            predicted_transformed = label_binarizer.transform(predicted)
            cv_acc = sklearn.metrics.accuracy_score(y_train, predicted_transformed)
            probabilities_d.update({act + '_' + str(nnum) + '_' + str(layer) + 'layers':probabilities})
            cv_scores.update({act + '_' + str(nnum) + '_' + str(layer) + 'layers':cv_acc})
            print('number of layers = %d; number of nodes = %d: %.3f' % (layer, nnum, cv_acc))
'''
with open('/mnt/isilon/diskin_lab/GMKF/tumor_data/deepsvr/innerlayers_cv_acc.txt','w') as outfile:
    outfile.write(str(cv_scores))
with open('/mnt/isilon/diskin_lab/GMKF/tumor_data/deepsvr/innerlayers_cv_probs.txt','w') as outfile:
    outfile.write(str(probabilities_d))
#np.savetxt('C:/Users/LEEL7/Documents/gmkf/deepsvr/test.csv',x_train,delimiter=',')
np.savetxt('/mnt/isilon/diskin_lab/GMKF/tumor_data/deepsvr/innerlayers_cv_xtrain.txt',x_train,delimiter=',')
np.savetxt('/mnt/isilon/diskin_lab/GMKF/tumor_data/deepsvr/innerlayers_cv_xtest.csv',x_test,delimiter=',')
np.savetxt('/mnt/isilon/diskin_lab/GMKF/tumor_data/deepsvr/innerlayers_cv_ytrain.csv',y_train,delimiter=',')
np.savetxt('/mnt/isilon/diskin_lab/GMKF/tumor_data/deepsvr/innerlayers_cv_ytest.csv',y_test,delimiter=',')

'''
#%%train??
#probabilities = cross_val_predict(nnEstimator, x_train, y_train, cv=kfold, method='predict_proba')
probabilities = probabilities_d['relu_50_3layers']

#%%classify test set
predicted = np.array([list(a).index(max(list(a))) for a in list(probabilities)])
#score the predictions
label_binarizer = preprocessing.LabelBinarizer()
label_binarizer.fit(range(max(predicted)+1))
predicted_transformed = label_binarizer.transform(predicted)


#%%
print('Cross validation accuracy:')
print('\t', sklearn.metrics.accuracy_score(y_train, predicted_transformed))
print('\nCross validation classification report\n')
print(sklearn.metrics.classification_report(y_train, predicted_transformed))

#%%

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


#%% Explore best parameters
#%% Determine the best network structure
'''
def evaluate_first_layer_build_cv_model(n_input,n_classes,nnodes,activationf):
    # define model
    cvmodel = keras.models.Sequential()
    cvmodel.add(Dense(nnodes, input_dim=n_input, activation=activationf, kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.001)))
    cvmodel.add(Dense(n_classes, activation='softmax'))
    cvmodel.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    return(cvmodel)
'''


def evaluate_first_layer(nnodes,activationf,trainx,trainy,testx,testy):
    
    # configure the model based on the data
    n_input, n_classes = trainx.shape[1], testy.shape[1]
    # define model
    model = keras.models.Sequential()
    model.add(Dense(nnodes, input_dim=n_input, activation=activationf, kernel_initializer='he_uniform',kernel_regularizer=keras.regularizers.l2(0.001)))
    model.add(Dense(n_classes, activation='softmax'))
    # compile model
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    print('performing cross-validation')
    # perform 10-fold CV
    cvscores = sklearn.model_selection.cross_val_score(KerasClassifier(build_fn=evaluate_first_layer_build_cv_model(n_input,n_classes,nnodes,activationf),epochs=1,
                                                                       batch_size=1500,verbose=1),x_train,y_train,cv=kfold)
    print('fitting model on training set')
    # fit model on train set
    fitmodel = model.fit(trainx, trainy, epochs=1, verbose=1)
    
    #fitmodel = model.fit(trainx, trainy, epochs=5, verbose=1)
    # evaluate model on test set
    _, test_acc = model.evaluate(testx, testy, verbose=1)
    return fitmodel, test_acc, cvscores

nnodes, activations = [50,75,100,125,150], ['relu','sigmoid','tanh']
#nnodes, activations = [125,150,200,250], ['relu']
histories = {}
test_results = {}
cv_scores = {}

for act in activations:
    print('activation function: {}'.format(act))
    for nnum in nnodes:
        # evaluate model with given params
        fitmodel, result, cvSC = evaluate_first_layer(nnum,act,x_train,y_train,x_test,y_test)
        histories.update({act + '_' + str(nnum):fitmodel.history})
        test_results.update({act + '_' + str(nnum):result})
        cv_scores.update({act + '_' + str(nnum):cvSC})
        print('number of nodes = %d: %.3f' % (nnum, result))

#%% plot results
linestyles = []
for key,item in histories.items():
    #print(key)
    plt.plot(item['loss'],label=key)
plt.legend()
plt.show()x
# show the plot
linestyles =  [(0, ()),
    (0, (1, 10)),
    (0, (1, 5)),
    (0, (1, 1)),
    (0, (5, 10)),
    (0, (5, 5)),
    (0, (5, 1)),
    (0, (3, 10, 1, 10)),
    (0, (3, 5, 1, 5)),
    (0, (3, 1, 1, 1)),
    (0, (3, 10, 1, 10, 1, 10)),
    (0, (3, 5, 1, 5, 1, 5)),
    (0, (3, 1, 1, 1, 1, 1))]
for i,(key,item) in enumerate(histories.items()):
    plt.plot(item['acc'],label=key,linestyle=linestyles[i%13])
plt.legend()
plt.show()

test_results
'''
    
    
    
