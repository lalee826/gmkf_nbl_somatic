#%% Project description
'''
- This script imports mutational signature decomposition results and clusters samples by their signature profiles using either a Gaussian mixture model or Kmeans
'''

import pandas as pd
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans
from sklearn import metrics
import numpy as np
from matplotlib import pyplot as plt

workdir = '/rocker-build/gmkf_nbl_somatic/'

#%% import data
d = pd.read_csv(workdir + 'Data/cosmicv3_sig_contribution.tsv',sep='\t',header=0,index_col=0)

##remove signatures for which max contribution is less than 5% for any sample
cnames = d.columns.values
toRemove = np.array([]).astype(str)
toRemove = np.concatenate((toRemove,(cnames[d.max(axis=0) <= 0.05]).astype(str)))
##remove signatures for which average contribution is less than 0.5%
toRemove = np.concatenate((toRemove,(cnames[d.mean(axis=0) < 0.005]).astype(str)))
##remove signatures for which max contribution is less than 10% for any sample, and over 90% of samples have 0% contribution
toRemove = np.concatenate((toRemove,(cnames[np.logical_and(d.max(axis=0) <= 0.1,d.apply(axis=0,func=lambda x:((np.sum(x==0))/348) >= 0.90))]).astype(str)))
#remove duplicates
toRemove = np.array(list(set(toRemove.tolist())))
len(toRemove)
keepCol = cnames[np.array([x not in toRemove for x in cnames])]
len(keepCol)
d = d[keepCol]
d.shape

#%% run kmeans and make distortion plot (testing up to 20 clusters)
kmod = [KMeans(n_clusters=k).fit(d) for k in np.arange(1,20)]
dst = [kmod[i].inertia_ for i in np.arange(len(kmod))]

#%% test clusters with a GMM and evaluate BIC and AIC (testing 20 iterations each for 2-15 clusters)
def SelBest(arr,x):
    '''
    returns the set of X configurations with shorter distance
    '''
    dx=np.argsort(arr)[:x]
    return arr[dx]
max_clusters = 15
nc = np.arange(2,max_clusters)
bics=[]
bics_err=[]
aics=[]
aics_err=[]
its=20
for n in nc:
    print(n)
    tmp_bic=[]
    tmp_aic=[]
    for x in range(its):
        gmmMod=GaussianMixture(n,n_init=2,covariance_type='full',random_state=0).fit(d) 
        tmp_bic.append(gmmMod.bic(d))
        tmp_aic.append(gmmMod.aic(d))
    valb=np.mean(SelBest(np.array(tmp_bic), int(its/5)))
    errb=np.std(tmp_bic)
    bics.append(valb)
    bics_err.append(errb)
    vala=np.mean(SelBest(np.array(tmp_aic), int(its/5)))
    erra=np.std(tmp_aic)
    aics.append(vala)
    aics_err.append(erra)

#%% plot bic and aic
plt.errorbar(nc, np.gradient(bics), yerr=bics_err, label='BIC')
plt.title("Gradient of BIC Scores", fontsize=20)
plt.xticks(nc)
plt.xlabel("N. of clusters")
plt.ylabel("grad(BIC)")
plt.legend()

#%%aic
plt.errorbar(nc, np.gradient(aics), yerr=aics_err, label='AIC')
plt.title("Gradient of AIC Scores", fontsize=20)
plt.xticks(nc)
plt.xlabel
("N. of clusters")
plt.ylabel("grad(AIC)")
plt.legend()

#%%GMM with silohuette score; same parameters as previous
n_clusters = np.arange(2,15)
sils = []
sils_err = []
its= 20
for n in n_clusters:
    print(n)
    tmp_sil = []
    #for each n_cluster perform 20 iterations
    for x in np.arange(its):
        gmmMod=GaussianMixture(n,n_init=2,covariance_type='full',random_state=0).fit(d) 
        labels=gmmMod.predict(d) #classify data based on fit model
        sil=metrics.silhouette_score(d,labels,metric='euclidean')#get silhouette score of model fit
        tmp_sil.append(sil)#append for iteration
    val=np.mean(SelBest(np.array(tmp_sil),int(its/5)))
    err=np.std(tmp_sil)
    sils.append(val)
    sils_err.append(err)

#%%plot
plt.errorbar(n_clusters, sils, yerr=sils_err)
plt.title("Silhouette Scores", fontsize=20)
plt.xticks(n_clusters)
plt.xlabel("N. of clusters")
plt.ylabel("Score")

#%%the kmeans distortion plot as well as the BIC/AIC/silhouette scores from GMM models suggest 8 final clusters be implemented
finCluster=8
mmMod=GaussianMixture(finCluster,n_init=5,covariance_type='full',random_state=0).fit(d) 
labels=mmMod.predict(d) #classify data based on fit model
ldf = pd.DataFrame({'sample_id':d.index.values,'cluster':labels})
ldf['cluster'] = [x+1 for x in ldf.cluster.tolist()]
#these cluster assignments are used in implementation of figure 4
#a results file is already provided in the Data folder