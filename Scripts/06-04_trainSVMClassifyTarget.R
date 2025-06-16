#%% Project description
'''
- This script performs SVM classification of TARGET gsva profiles using the GMKF assigned classifications as training.
- 50 iterations are run and the mode group is assigned to each TARGET sample
- Visualization performed with either umap or tsne
- Sample order for TARGET is exported to be plotted in R
'''

import pandas as pd
import numpy as np
from sklearn import svm
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from statistics import mode
import matplotlib.pyplot as plt
import plotly.express as px
import plotly
plotly.io.renderers.default='browser'
workdir = '/rocker-build/gmkf_nbl_somatic/'

#%% import data and prepare data frames for input
gFilePath = workdir + 'Data/gmkf_gsva_for_classification.tsv'
sFilePath = workdir + 'Data/target_gsva_for_classification.tsv'
gData = pd.read_csv(gFilePath,sep='\t',header=0).set_index('case_id')
sData = pd.read_csv(sFilePath,sep='\t',header=0).set_index('case_id')
gClass = np.array(gData['class'])
gData = gData.drop(['class'],axis=1)

#import clinical data
cdgFilePath = workdir + 'Data/cognbl_gmkf_clindata.csv'
sdgFilePath = workdir + 'Data/TARGET_harmonized_2018-03-31.csv'
cdg = pd.read_csv(cdgFilePath,sep=',',header=0)
sdg = pd.read_csv(sdgFilePath,sep=',',header=0)

#%% scale each feature to standard normal
stdScaler = StandardScaler()
gDataScaled = stdScaler.fit_transform(gData)
sDataScaled = stdScaler.fit_transform(sData)

#%% fit gkmf data into SVM model
svmModel = svm.SVC( C=1.0, kernel='rbf', degree=3, gamma='scale', coef0=0.0, shrinking=True, probability=False, tol=0.001, cache_size=500, 
                   class_weight=None, verbose=False, max_iter=-1, decision_function_shape='ovr', random_state=None)
svmModel.fit(gDataScaled,gClass)

#%% predict classification of TARGET data
preds = svmModel.predict(sDataScaled)

#%% run multiple iterations with random regularization and take mode of predictions for each sample
iterations = 150
kernel_type = 'rbf'

predictArray = np.empty((0,sDataScaled.shape[0]),int)

for i in np.arange(1,iterations+1):
    #use variable regularization
    regVal = int(np.random.uniform(1,10,size=1)[0])
    svmModel = svm.SVC( C=regVal, kernel=kernel_type, degree=3, gamma='scale', coef0=0.0, shrinking=True, probability=False, tol=0.001, cache_size=500, 
                   class_weight=None, verbose=False, max_iter=-1, decision_function_shape='ovr', random_state=None)
    svmModel.fit(gDataScaled,gClass)
    predsi = svmModel.predict(sDataScaled)
    predictArray = np.append(predictArray,np.array([predsi]),axis=0)

finalPreds = np.apply_along_axis(mode,axis=0,arr=predictArray)

#%% visualize using t-SNE or UMAP
#plot in 2d
method = 'raw'
if method == 'scaled':
    plotData = pd.concat([pd.DataFrame(gDataScaled),pd.DataFrame(sDataScaled)])
else:
    plotData = pd.concat([gData,sData])
#plotData = pd.concat([pd.DataFrame(gDataScaled),pd.DataFrame(sDataScaled)])
#testing shows ideal parameters are mindist 0.1 and spread 2-4
minDist = 0.01
n_neighbors = 10
spread=5
#umapModel = UMAP(n_neighbors=n_neighbors,n_components=2,n_epochs=500,metric='euclidean',min_dist=minDist,spread=spread)
tsneModel = TSNE(n_components=2, perplexity=20, early_exaggeration=12.0, learning_rate=200.0, n_iter=1000, n_iter_without_progress=300, 
                 min_grad_norm=1e-07, metric='euclidean', init='random', verbose=0)
plotCords = tsneModel.fit_transform(plotData)

#plot in 2d using plotly
ng = len(gClass)
ns = len(finalPreds)
nt = plotData.shape[0]
plotCordsDF = pd.DataFrame(plotCords).rename({0:'X',1:'Y'},axis=1)
plotCordsDF['class'] = [str(x) for x in gClass] + [str(x) for x in finalPreds]
plotCordsDF['size'] = [50]*ng + [12]*ns
plotCordsDF['symbol'] = ['GMKF']*ng + ['TARGET']*ns
plotCordsDF['opacity'] = [0.7]*nt
getSurv = np.vectorize(lambda x: cdg.loc[cdg.usi==x,'vital_status'].tolist()[0] if x in gData.index.values else sdg.loc[sdg.usi==x,'Vital Status'].tolist()[0])
plotCordsDF['survival'] = getSurv(np.concatenate([gData.index.values,sData.index.values]))
plotCordsDF['group'] = [plotCordsDF.loc[i,'symbol'] + '-' + plotCordsDF.loc[i,'survival'] + '-' + plotCordsDF.loc[i,'class'] for i in np.arange(0,plotCordsDF.shape[0])]

##%create a separate trace for each factor level
#create a color for each trace
color_map = {
    1:'yellow',
    2:'darkorange',
    3:'darkgreen',
    4:'darkblue',
    5:'lawngreen',
    6:'purple',
    7:'aquamarine',
    8:'darksalmon',
    9:'slategray',
    10:'mediumorchid',
    11:'maroon',
}
#set parameters for each trace
factorMap = {
    'GMKF-alive-1':{'outline':'black','size':35,'shape':'circle','fill':color_map[1],'opacity':0.8},
    'GMKF-alive-2':{'outline':'black','size':35,'shape':'circle','fill':color_map[2],'opacity':0.8},
    'GMKF-alive-3':{'outline':'black','size':35,'shape':'circle','fill':color_map[3],'opacity':0.8},
    'GMKF-alive-4':{'outline':'black','size':35,'shape':'circle','fill':color_map[4],'opacity':0.8},
    'GMKF-alive-5':{'outline':'black','size':35,'shape':'circle','fill':color_map[5],'opacity':0.8},
    'GMKF-alive-6':{'outline':'black','size':35,'shape':'circle','fill':color_map[6],'opacity':0.8},
    'GMKF-alive-7':{'outline':'black','size':35,'shape':'circle','fill':color_map[7],'opacity':0.8},
    'GMKF-alive-8':{'outline':'black','size':35,'shape':'circle','fill':color_map[8],'opacity':0.8},
    'GMKF-alive-9':{'outline':'black','size':35,'shape':'circle','fill':color_map[9],'opacity':0.8},
    #'GMKF-alive-10':{'outline':'black','size':35,'shape':'circle','fill':color_map[10],'opacity':0.8},
    #'GMKF-alive-11':{'outline':'black','size':35,'shape':'circle','fill':color_map[11],'opacity':0.8},
    #'GMKF-dead-1':{'outline':'red','size':35,'shape':'circle-dot','fill':color_map[1],'opacity':0.8},
    'GMKF-dead-2':{'outline':'red','size':35,'shape':'circle-dot','fill':color_map[2],'opacity':0.8},
    'GMKF-dead-3':{'outline':'red','size':35,'shape':'circle-dot','fill':color_map[3],'opacity':0.8},
    'GMKF-dead-4':{'outline':'red','size':35,'shape':'circle-dot','fill':color_map[4],'opacity':0.8},
    #'GMKF-dead-5':{'outline':'red','size':35,'shape':'circle-dot','fill':color_map[5],'opacity':0.8},
    #'GMKF-dead-6':{'outline':'red','size':35,'shape':'circle-dot','fill':color_map[6],'opacity':0.8},
    'GMKF-dead-7':{'outline':'red','size':35,'shape':'circle-dot','fill':color_map[7],'opacity':0.8},
    #'GMKF-dead-8':{'outline':'red','size':35,'shape':'circle-dot','fill':color_map[8],'opacity':0.8},
    'GMKF-dead-9':{'outline':'red','size':35,'shape':'circle-dot','fill':color_map[9],'opacity':0.8},
    #'GMKF-dead-10':{'outline':'red','size':35,'shape':'circle-dot','fill':color_map[10],'opacity':0.8},
    #'GMKF-dead-11':{'outline':'red','size':35,'shape':'circle-dot','fill':color_map[11],'opacity':0.8},
    'TARGET-Alive-1':{'outline':'black','size':20,'shape':'diamond','fill':color_map[1],'opacity':0.65},
    'TARGET-Alive-2':{'outline':'black','size':20,'shape':'diamond','fill':color_map[2],'opacity':0.65},
    'TARGET-Alive-3':{'outline':'black','size':20,'shape':'diamond','fill':color_map[3],'opacity':0.65},
    'TARGET-Alive-4':{'outline':'black','size':20,'shape':'diamond','fill':color_map[4],'opacity':0.65},
    'TARGET-Alive-5':{'outline':'black','size':20,'shape':'diamond','fill':color_map[5],'opacity':0.65},
    'TARGET-Alive-6':{'outline':'black','size':20,'shape':'diamond','fill':color_map[6],'opacity':0.65},
    'TARGET-Alive-7':{'outline':'black','size':20,'shape':'diamond','fill':color_map[7],'opacity':0.65},
    #'TARGET-Alive-8':{'outline':'black','size':20,'shape':'diamond','fill':color_map[8],'opacity':0.65},
    'TARGET-Alive-9':{'outline':'black','size':20,'shape':'diamond','fill':color_map[9],'opacity':0.65},
    #'SEQC-0-10':{'outline':'black','size':20,'shape':'diamond','fill':color_map[10],'opacity':0.65},
    #'SEQC-0-11':{'outline':'black','size':20,'shape':'diamond','fill':color_map[11],'opacity':0.65},
    'TARGET-Dead-1':{'outline':'red','size':20,'shape':'diamond-dot','fill':color_map[1],'opacity':0.65},
    'TARGET-Dead-2':{'outline':'red','size':20,'shape':'diamond-dot','fill':color_map[2],'opacity':0.65},
    'TARGET-Dead-3':{'outline':'red','size':20,'shape':'diamond-dot','fill':color_map[3],'opacity':0.65},
    'TARGET-Dead-4':{'outline':'red','size':20,'shape':'diamond-dot','fill':color_map[4],'opacity':0.65},
    'TARGET-Dead-5':{'outline':'red','size':20,'shape':'diamond-dot','fill':color_map[5],'opacity':0.65},
    'TARGET-Dead-6':{'outline':'red','size':20,'shape':'diamond-dot','fill':color_map[6],'opacity':0.65},
    'TARGET-Dead-7':{'outline':'red','size':20,'shape':'diamond-dot','fill':color_map[7],'opacity':0.65},
    #'TARGET-Dead-8':{'outline':'red','size':20,'shape':'diamond-dot','fill':color_map[8],'opacity':0.65},
    'TARGET-Dead-9':{'outline':'red','size':20,'shape':'diamond-dot','fill':color_map[9],'opacity':0.65},
    #'SEQC-1-10':{'outline':'red','size':20,'shape':'diamond-dot','fill':color_map[10],'opacity':0.65},
    #'SEQC-1-11':{'outline':'red','size':20,'shape':'diamond-dot','fill':color_map[11],'opacity':0.65}
}

#create plotly traces
traces = []
for factor, setting in factorMap.items():
    print(factor)
    print(setting)
    trace = px.scatter(plotCordsDF.loc[plotCordsDF.group==factor],x='X',y='Y',
                       color='class',#color_discrete_map=color_map,
                       #size='size',#size_sequence=[szMap[factor]],
                       symbol='symbol',symbol_sequence=[factorMap[factor]['shape']],
                       opacity=factorMap[factor]['opacity']).update_traces(marker=dict(line=dict(width=3, color=factorMap[factor]['outline']),
                                                              size=factorMap[factor]['size'],
                                                              color=factorMap[factor]['fill'])
                                                              )
    for i in np.arange(0,len(trace.data)):
        print(i)
        traces.append(trace.data[i])

#initialize figure
fig = px.scatter(title='mindist:' + str(minDist) + ';\nspread:' + str(spread) + ';\nneighbors:' + str(n_neighbors) + 
                 ';\ndata:' + method + ';\nsvm:' + kernel_type,)

#add each trace
for trace in traces:
    fig.add_trace(trace)
fig.show()

#%% check survival between data sets for each group
survDF = pd.DataFrame({'group':np.arange(1,9),
                       'GMKF_alive':[plotCordsDF.loc[np.logical_and(np.logical_and(plotCordsDF.symbol=='GMKF',plotCordsDF['class']==str(x)),
                                                                    plotCordsDF.survival=='alive')].shape[0] for x in np.arange(1,9)],
                       'GMKF_dead':[plotCordsDF.loc[np.logical_and(np.logical_and(plotCordsDF.symbol=='GMKF',plotCordsDF['class']==str(x)),
                                                                    plotCordsDF.survival=='dead')].shape[0] for x in np.arange(1,9)],
                       'TARGET_alive':[plotCordsDF.loc[np.logical_and(np.logical_and(plotCordsDF.symbol=='TARGET',plotCordsDF['class']==str(x)),
                                                                    plotCordsDF.survival=='Alive')].shape[0] for x in np.arange(1,9)],
                       'TARGET_dead':[plotCordsDF.loc[np.logical_and(np.logical_and(plotCordsDF.symbol=='TARGET',plotCordsDF['class']==str(x)),
                                                                    plotCordsDF.survival=='Dead')].shape[0] for x in np.arange(1,9)]                      
                      })

#%% check by group
checkG2 = survDF.loc[1]

#%% sample order for TARGET
sosc = plotCordsDF.loc[plotCordsDF.symbol=='TARGET']['class'].reset_index(drop=True).rename(
    index=dict(tuple(zip(np.arange(0,sdg.shape[0]),sData.index.values.tolist())))).reset_index()
sosc.groupby('class').count()