from sklearn.datasets.samples_generator import make_blobs
import os
from sklearn.cluster import KMeans
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from copy import deepcopy
import time
pd.options.display.max_rows = 500




def make_array(x):
    arr = []
    for i in range(1,x+1):
        arr.append(i)
    return arr


def test_clust():
     os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
     phipsidf = pd.read_hdf("phipsi_xray_chains_hetatmlisted.h5",)
     clustlist = [[-69.611298,160.439148,67.632287,28.367751,-101.183084,-12.240481,86.874515,11.805842,-111.793021,20.48053,63.360547,39.072422,-58.639519,-50.749508],
                  [-91.549515,146.792979,50.120174,29.000047,-119.388986,41.483409,74.403471,11.448555,-90.089563,-7.485303,77.770843,40.993155,-57.168942,-36.872158],
                  [-94.809597,	143.357633,63.971899,	33.944542,-122.715199	,7.485939,91.586699	,14.426918,-105.533392	,-9.983716,70.332541	,19.947048,-53.669792	,-44.714479],
                  [-73.925067,	155.701969,49.182985,	40.227576,-124.332632,	12.530059,82.147041,	9.776688,-100.668852	,15.728854,63.512251,	39.026531,-60.964152	,-35.683852],
                  [-83.568937,	130.098797,68.287428	,21.33896,-88.304481	,-3.263119,87.716364	,16.252141,-105.853914	,16.741628,71.005686	,53.915787,-51.255142	,-44.249912],
                  [-82.399156,	160.907402,84.038429	,12.187862,-102.244899	,4.391591,55.080186	,43.018486,-135.32258	,9.678417,59.326935	,32.669506,-62.596224,	-38.157732]]
     testdf = pd.DataFrame(clustlist, columns=[1,2,3,4,5,6,7,8,9,10,11,12,13,14])
     #print(testdf)
     kmeans = KMeans(n_clusters = 3)
     kmeans.fit(testdf)
     labels = kmeans.predict(testdf)
     centroids = kmeans.cluster_centers_
     testdf[15] = labels
     print(testdf)

def get_phipsi(line_number,namechainres,length):
    name = (namechainres)[0:4]
    chain = (namechainres)[10:11]
    res= (namechainres)[15:]
    line_number = line_number-1
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/phipsi_xray_chains")
    filename = "pdb"+name+"_chain"+chain+"_biopython.tsv"
    angledf = pd.read_csv(filename, sep='\t',usecols=[0,1,2], names=['name_chain_res','phi','psi'])
    angle_list = []
    for index, row in angledf.iloc[line_number:line_number+length].iterrows():
        angle_list.append(row['phi'])
        angle_list.append(row['psi'])
    #angle_list = np.array(angle_list)
    return angle_list

def get_phipsi_trig(line_number,namechainres,length):
    name = (namechainres)[0:4]
    chain = (namechainres)[10:11]
    res= (namechainres)[15:]
    line_number = line_number-1
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/phipsi_xray_chains")
    filename = "pdb"+name+"_chain"+chain+"_biopython.tsv"
    angledf = pd.read_csv(filename, sep='\t',usecols=[0,1,2], names=['name_chain_res','phi','psi'])
    angle_list = []
    for index, row in angledf.iloc[line_number:line_number+length].iterrows():
        angle_list.append(np.cos(np.deg2rad(row['phi'])))
        angle_list.append(np.sin(np.deg2rad(row['phi'])))
        angle_list.append(np.cos(np.deg2rad(row['psi'])))
        angle_list.append(np.sin(np.deg2rad(row['psi'])))
    #angle_list = np.array(angle_list)
    return angle_list


def get_phipsi_get_pos(line_number,namechainres,length):
    name = (namechainres)[0:4]
    chain = (namechainres)[10:11]
    res= (namechainres)[15:]
    line_number = line_number-1
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/phipsi_xray_chains")
    filename = "pdb"+name+"_chain"+chain+"_biopython.tsv"
    angledf = pd.read_csv(filename, sep='\t',usecols=[0,1,2], names=['name_chain_res','phi','psi'])
    angle_list = []
    for index, row in angledf.iloc[line_number:line_number+length].iterrows():
        angle_list.append(row['phi'])
        angle_list.append(row['psi'])
    #angle_list = np.array(angle_list)
    if angle_list[0] < 0.0:
        angle_list.append(0)
    else:
        angle_list.append(1)
    return angle_list

def get_phipsi_get_pos_trig(line_number,namechainres,length):
    name = (namechainres)[0:4]
    chain = (namechainres)[10:11]
    res= (namechainres)[15:]
    line_number = line_number-1
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/phipsi_xray_chains")
    filename = "pdb"+name+"_chain"+chain+"_biopython.tsv"
    angledf = pd.read_csv(filename, sep='\t',usecols=[0,1,2], names=['name_chain_res','phi','psi'])
    angle_list = []
    pos_list = []
    for index, row in angledf.iloc[line_number:line_number+length].iterrows():
        pos_list.append(row['phi'])
        angle_list.append(np.cos(np.deg2rad(row['phi'])))
        angle_list.append(np.sin(np.deg2rad(row['phi'])))
        angle_list.append(np.cos(np.deg2rad(row['psi'])))
        angle_list.append(np.sin(np.deg2rad(row['psi'])))
    #angle_list = np.array(angle_list)
    if pos_list[0] < 0.0:
        angle_list.append(0)
    else:
        angle_list.append(1)
    return angle_list

def get_clustlist(length):
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    phipsidf = pd.read_csv("phi_xray_chains_sorted_curvature_hetatm_rev4.csv",usecols = ['line_number','name_chain_res','length','curvature','hetatm'])
    #phipsidf = phipsidf.sort_values(by ='length', ascending = False)
    dfL = phipsidf[phipsidf['length']==length]  
    #print(dfL)
    df_temp =  pd.DataFrame()
    #df_temp[make_array(14)]= zip(*df7[['line_number','name_chain_res','length']].map(get_phipsi_apply))
    df_temp[make_array(length*2)]= dfL.apply(lambda row: get_phipsi(row['line_number'], row['name_chain_res'], row['length']),axis = 1,result_type='expand')
    #print(df_temp)
    return dfL,df_temp

def get_clustlist_get_pos(length):
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    phipsidf = pd.read_csv("phi_xray_chains_sorted_curve_bfac_is_pos_cluster_conf_hetatm_rev4-20clust.csv",usecols = ['line_number','name_chain_res','length','curvature','hetatm'])
    #phipsidf = phipsidf.sort_values(by ='length', ascending = False)
    dfL = deepcopy(phipsidf[phipsidf['length']==length])
    df_data =  pd.DataFrame()
    df_data[make_array((length*2)+1)]= dfL.apply(lambda row: get_phipsi_get_pos(row['line_number'], row['name_chain_res'], row['length']),axis = 1,result_type='expand')
    dfL['is_pos'] = df_data[(length*2)+1]
    dfL = dfL.sort_values(by = 'is_pos',ascending=False)
    df_data = df_data.sort_values(by = ((length*2)+1), ascending = False)

    return dfL,df_data

def get_clustlist_get_pos_trig(length):
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    phipsidf = pd.read_csv("phi_xray_chains_sorted_curve_bfac_is_pos_cluster_conf_hetatm_rev4-20clust.csv",usecols = ['line_number','name_chain_res','length','bfac_within_stddev','curvature','hetatm'])
    #phipsidf = phipsidf.sort_values(by ='length', ascending = False)
    dfL = deepcopy(phipsidf[phipsidf['length']==length])
    df_data =  pd.DataFrame()
    df_data[make_array((length*4)+1)]= dfL.apply(lambda row: get_phipsi_get_pos_trig(row['line_number'], row['name_chain_res'], row['length']),axis = 1,result_type='expand')
    dfL['is_pos'] = df_data[(length*4)+1]
    dfL = dfL.sort_values(by = 'is_pos',ascending=False)
    df_data = df_data.sort_values(by = ((length*4)+1), ascending = False)

    return dfL,df_data

def get_clustlist_trig(length):
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    phipsidf = pd.read_csv("phi_xray_chains_sorted_curve_bfac_is_pos_cluster_conf_hetatmmulti_mse_edo_unx_25clust_TRIG_rev4-20clust.csv",usecols = ['line_number','name_chain_res','length','bfac_within_stddev','curvature','is_pos','cluster','conf_seq','hetatm','hetatm_simple'])
    #phipsidf = phipsidf.sort_values(by ='length', ascending = False)
    dfL = deepcopy(phipsidf[phipsidf['length']==length])
    df_data =  pd.DataFrame()
    df_data[make_array((length*4))]= dfL.apply(lambda row: get_phipsi_trig(row['line_number'], row['name_chain_res'], row['length']),axis = 1,result_type='expand')
    dfL['is_pos'] = df_data[(length*4)]
    dfL = dfL.sort_values(by = 'is_pos',ascending=False)
    df_data = df_data.sort_values(by = ((length*4)+1), ascending = False)

    return dfL,df_data


def big_clust(clustdf,length,k):
     df_temp = clustdf[1]
     kmeans = KMeans(n_clusters = k)
     kmeans.fit(df_temp)
     labels = kmeans.predict(df_temp)
     centroids = kmeans.cluster_centers_
     df_temp['cluster'] = labels
     print(df_temp)
     dfL = clustdf[0]
     dfL['cluster'] = labels
     dfL = dfL.sort_values(by = ['is_pos', 'cluster'],ascending=False)
     return dfL





def big_clust_get_pos(clustdf,length,k1,k0):
    df_data = clustdf[1]
    dfL = clustdf[0]
    df_data1 = deepcopy(df_data[df_data[((length*2)+1)]==1])
    dfL1 = deepcopy(dfL[dfL['is_pos']==1])
    df_data0 = deepcopy(df_data[df_data[((length*2)+1)]==0])
    dfL0 = deepcopy(dfL[dfL['is_pos']==0])
    kmeans1 = KMeans(n_clusters = k1)
    kmeans1.fit(df_data1)
    labels1 = kmeans1.predict(df_data1)
    centroids = kmeans1.cluster_centers_
    df_data1['cluster'] = labels1
    #print(df_data1)
    kmeans0 = KMeans(n_clusters = k0)
    kmeans0.fit(df_data0)
    labels1 = kmeans0.predict(df_data0)
    centroids = kmeans0.cluster_centers_
    df_data0['cluster'] = labels1
    #print(df_data0)    
    df_data_total = pd.concat([df_data1,df_data0])
    print(df_data_total)
    dfL['cluster'] = df_data_total['cluster']
    dfL = dfL.sort_values(by = ['is_pos', 'cluster'],ascending=False)
    return dfL

def big_clust_get_pos_trig(clustdf,length,k1,k0):
    df_data = clustdf[1]
    dfL = clustdf[0]
    df_data1 = deepcopy(df_data[df_data[((length*4)+1)]==1])
    dfL1 = deepcopy(dfL[dfL['is_pos']==1])
    df_data0 = deepcopy(df_data[df_data[((length*4)+1)]==0])
    dfL0 = deepcopy(dfL[dfL['is_pos']==0])
    kmeans1 = KMeans(n_clusters = k1)
    kmeans1.fit(df_data1)
    labels1 = kmeans1.predict(df_data1)
    centroids = kmeans1.cluster_centers_
    df_data1['cluster'] = labels1
    #print(df_data1)
    kmeans0 = KMeans(n_clusters = k0)
    kmeans0.fit(df_data0)
    labels1 = kmeans0.predict(df_data0)
    centroids = kmeans0.cluster_centers_
    df_data0['cluster'] = labels1
    #print(df_data0)    
    df_data_total = pd.concat([df_data1,df_data0])
    print(df_data_total)
    dfL['cluster'] = df_data_total['cluster']
    dfL = dfL.sort_values(by = ['is_pos', 'cluster'],ascending=False)
    return dfL


def elbow(clustdf):
    ssd = []
    krange = range(1,5)
    for k in krange:
        kmeans = KMeans(n_clusters = k)
        kmeans.fit(clustdf[1])
        ssd.append(kmeans.inertia_)
    plt.style.use('seaborn')
    plt.plot(krange, ssd, 'bx-')
    plt.xlabel('k')
    plt.ylabel('Sum of squared distances')
    plt.title('Elbow Method For Optimal k')
    plt.show()  

def elbow_is_pos(clustdf,length):
    df_data = clustdf[1]
    dfL = clustdf[0]
    df_data1 = deepcopy(df_data[df_data[((length*2)+1)]==1])
    dfL1 = deepcopy(dfL[dfL['is_pos']==1])
    print(df_data1,dfL1)
    df_data0 = deepcopy(df_data[df_data[((length*2)+1)]==0])
    dfL0 = deepcopy(dfL[dfL['is_pos']==0])
    print(df_data0,dfL0)
    ssd = []
    krange = range(1,5)
    for k in krange:
        kmeans = KMeans(n_clusters = k)
        kmeans.fit(df_data1)
        ssd.append(kmeans.inertia_)
    plt.style.use('seaborn')
    plt.plot(krange, ssd, 'bx-')
    plt.xlabel('k')
    plt.ylabel('Sum of squared distances')
    plt.title('Elbow Method For Optimal k')
    plt.show()  



def elbow_is_pos_trig(clustdf,length):
    df_data = clustdf[1]
    dfL = clustdf[0]
    df_data1 = deepcopy(df_data[df_data[((length*4)+1)]==1])
    dfL1 = deepcopy(dfL[dfL['is_pos']==1])
    print(df_data1,dfL1)
    df_data0 = deepcopy(df_data[df_data[((length*4)+1)]==0])
    dfL0 = deepcopy(dfL[dfL['is_pos']==0])
    print(df_data0,dfL0)
    ssd = []
    krange = range(1,7)
    for k in krange:
        kmeans = KMeans(n_clusters = k)
        kmeans.fit(df_data0)
        ssd.append(kmeans.inertia_)
    plt.style.use('seaborn')
    plt.plot(krange, ssd, 'bx-')
    plt.xlabel('k')
    plt.ylabel('Sum of squared distances')
    plt.title('Elbow Method For Optimal k')
    plt.show()  
        
#df4 = big_clust(get_clustlist(4),4,2)
#df5 = big_clust(get_clustlist(5),5,3)
#df6 = big_clust(get_clustlist(6),6,2)
#df7 = big_clust(get_clustlist(7),7,2)
#dflist = [df7,df6,df5,df4]
#bigdf = pd.concat(dflist)
#os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
#bigdf.to_csv('phipsi_xray_chains_sorted_curvature_clustered_rev2.csv')
#print(bigdf)

#df4 = get_clustlist_get_pos(4)
#df5 = get_clustlist_get_pos(5)
#df6 = get_clustlist_get_pos(6)
#df7 = get_clustlist_get_pos(7)
#df8 = get_clustlist_get_pos(8)
#df9 = get_clustlist_get_pos(9)
#df11 = get_clustlist_get_pos(11)
#dflist = [df11,df9,df8,df7,df6,df5,df4]
#bigdf = pd.concat(dflist)
#os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
#bigdf.to_csv('phi_xray_chains_sorted_curvature_is_pos_cluster_hetatm_rev4.csv')
#print(bigdf)

#elbow(get_clustlist_get_pos(9))
#elbow_is_pos(get_clustlist_get_pos(5),5)

#print(get_clustlist_get_pos_trig(8)[0])

def main_cluster_runner():
    df4 = big_clust_get_pos_trig(get_clustlist_get_pos_trig(4),4,5,5)
    df5 = big_clust_get_pos_trig(get_clustlist_get_pos_trig(5),5,2,25)
    df6 = big_clust_get_pos_trig(get_clustlist_get_pos_trig(6),6,6,4)
    df7 = big_clust(get_clustlist_get_pos_trig(7),7,3)
    df8 = big_clust_get_pos_trig(get_clustlist_get_pos_trig(8),8,2,2)
    df9 = big_clust(get_clustlist_get_pos_trig(9),9,2)
    df11 = big_clust(get_clustlist_get_pos_trig(11),11,2)
    bigdf = pd.concat([df11,df9,df8,df7,df6,df5,df4])
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    bigdf.to_csv('phi_xray_chains_sorted_curve_bfac_is_pos_cluster_conf_hetatm_25clust_TRIG_rev4.csv')
    print(bigdf)



def get_subclustlist_get_pos_trig(length,is_pos,cluster):
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    phipsidf = pd.read_csv("phi_xray_chains_sorted_curve_bfac_is_pos_cluster_conf_hetatmmulti_mse_edo_unx_25clust_TRIG_rev4.csv",usecols = ['line_number','name_chain_res','length','bfac_within_stddev','curvature','is_pos','cluster','conf_seq','hetatm','hetatm_simple'])
    dfL = deepcopy(phipsidf[(phipsidf['length']==length) & (phipsidf['is_pos'] ==is_pos) & (phipsidf['cluster']==cluster)])
    df_data =  pd.DataFrame()
    df_data[make_array((length*4)+1)]= dfL.apply(lambda row: get_phipsi_get_pos_trig(row['line_number'], row['name_chain_res'], row['length']),axis = 1,result_type='expand')
    #dfL['is_pos'] = df_data[(length*4)+1] #dont need to modify dfL since is_pos is already there
    #dfL = dfL.sort_values(by = 'is_pos',ascending=False) #not gonna do anything since is_pos is only 0
    #df_data = df_data.sort_values(by = ((length*4)+1), ascending = False) #also not going to do anything for the same reason

    return dfL,df_data

def big_subclust(clustdf,length,k):
     df_temp = clustdf[1]
     kmeans = KMeans(n_clusters = k)
     kmeans.fit(df_temp)
     labels = kmeans.predict(df_temp)
     centroids = kmeans.cluster_centers_
     df_temp['subcluster'] = labels
     #print(df_temp)
     dfL = clustdf[0]
     dfL['subcluster'] = labels
     dfL = dfL.sort_values(by = 'subcluster',ascending=False)
     return dfL



def subcluster_runner():
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    phipsidf = pd.read_csv("phi_xray_chains_sorted_curve_bfac_is_pos_cluster_conf_hetatmmulti_mse_edo_unx_25clust_TRIG_rev4.csv",usecols = ['line_number','name_chain_res','length','bfac_within_stddev','curvature','is_pos','cluster','conf_seq','hetatm','hetatm_simple'])
    phipsidf['subcluster'] = 0
    df503 = big_subclust(get_subclustlist_get_pos_trig(5,0,3),5,3)
    df505 = big_subclust(get_subclustlist_get_pos_trig(5,0,5),5,2)
    df5023 =big_subclust(get_subclustlist_get_pos_trig(5,0,23),5,2)
    df5021 =big_subclust(get_subclustlist_get_pos_trig(5,0,21),5,2)
    #phipsidf = phipsidf[((phipsidf.length != 5) & (phipsidf.is_pos != 0) & (phipsidf.cluster != 3))] #cut out anything thats not 503,5,23,21
    #phipsidf = phipsidf[((phipsidf.length != 5) & (phipsidf.is_pos != 0) & (phipsidf.cluster != 5))]
    #phipsidf = phipsidf[((phipsidf.length != 5) & (phipsidf.is_pos != 0) & (phipsidf.cluster != 23))]
    #phipsidf = phipsidf[((phipsidf.length != 5) & (phipsidf.is_pos != 0) & (phipsidf.cluster != 21))]
    phipsidf.drop(phipsidf[(phipsidf['length']==5) & (phipsidf['is_pos'] == 0) & (phipsidf['cluster'] == 3)].index, inplace = True) #cut out 503 etc
    phipsidf.drop(phipsidf[(phipsidf['length']==5) & (phipsidf['is_pos'] == 0) & (phipsidf['cluster'] == 5)].index, inplace = True)
    phipsidf.drop(phipsidf[(phipsidf['length']==5) & (phipsidf['is_pos'] == 0) & (phipsidf['cluster'] == 23)].index, inplace = True)
    phipsidf.drop(phipsidf[(phipsidf['length']==5) & (phipsidf['is_pos'] == 0) & (phipsidf['cluster'] == 21)].index, inplace = True)

    print(phipsidf)
    #append new pieces then sort
    subclustlist = [df503,df505,df5023,df5021]
    subdf = pd.concat(subclustlist)
    phipsidf = pd.concat([subdf,phipsidf])
    phipsidf = phipsidf.sort_values(['length','is_pos','cluster','subcluster'],ascending = [False,False,False,False])
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    phipsidf.to_csv("phi_xray_chains_sorted_curve_bfac_is_pos_cluster_subcluster_conf_hetatmmulti_mse_edo_unx_25clust_TRIG_rev4.csv")
    #print(phipsidf)


     

    




start = time.time()
#good numbers of subclusters to use
#print(big_subclust(get_subclustlist_get_pos_trig(5,0,3),5,3))
#print(big_subclust(get_subclustlist_get_pos_trig(5,0,5),5,2))
#print(big_subclust(get_subclustlist_get_pos_trig(5,0,23),5,2))
#print(big_subclust(get_subclustlist_get_pos_trig(5,0,21),5,2))
#print(big_subclust(get_subclustlist_get_pos_trig(5,0,22),5,1))
 
#elbow_is_pos_trig(get_subclustlist_get_pos_trig(5,0,5),5)
subcluster_runner()
end = time.time()
print ('Run time is:'+str(end - start))


#df4 = big_clust_get_pos(get_clustlist_get_pos(4),4,3,3)
#df5 = big_clust_get_pos(get_clustlist_get_pos(5),5,2,4)
#df6 = big_clust_get_pos(get_clustlist_get_pos(6),6,2,2)
#df7 = big_clust(get_clustlist_get_pos(7),7,2)
#bigdf = pd.concat([df7,df6,df5,df4])
#os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
#bigdf.to_csv('phipsi_xray_chains_sorted_curvature_isPos_clustered_hetatm_rev2.csv')
#print(bigdf)


