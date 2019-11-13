import os
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from copy import deepcopy
import ast
import itertools as it


def split_hetstats(length,is_pos,cluster):
     os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
     phipsidf = pd.read_csv("phi_xray_chains_sorted_curve_bfac_is_pos_cluster_conf_hetatmmulti_mse_edo_unx_25clust_TRIG_rev4.csv",usecols = ['line_number','name_chain_res','length','bfac_within_stddev','curvature','is_pos','cluster','conf_seq','hetatm','hetatm_simple'])
     dfL = deepcopy(phipsidf[(phipsidf['length']==length) & (phipsidf['is_pos'] ==is_pos) & (phipsidf['cluster']==cluster)])
     df_blanks= deepcopy(dfL[(dfL['hetatm_simple']=='[]')])
     #print(df_blanks)
     percent_het = float(len(dfL)-len(df_blanks))/len(dfL)
     #print(percent_het)
     hetlist = []
     #hetlist.append(ast.literal_eval(dfL['hetatm simple']))
     for index,row in dfL.iterrows():
         hets = ast.literal_eval(row['hetatm_simple'])
         for het in hets:
            hetlist.append(het[0])
     hetdict = {}
     for het in hetlist:
         if het in hetdict:
             hetdict[het] +=1
         else:
            hetdict[het] = 1
     #print(hetdict)
     hetdict = [(hetdict[key], key) for key in hetdict]
     hetdict.sort()
     hetdict.reverse()
     top_hets = []
     added_data = [percent_het]
     for i in range(3):         
         added_data.append([hetdict[i][1],hetdict[i][0]])
     #print(top_het)
     return added_data

def combine_hetstats():
     df_export = pd.read_csv("phi_xray_chains_sorted_curve_bfac_is_pos_cluster_conf_hetatmstats_mse_25clust_trig_SEQS_rev4.csv",usecols = ['length','is_pos', 'cluster','# entries','avg_curvature','consensus'])
     added_data =[]
     for i in range(25):
         added_data.append(split_hetstats(5,0,i))
     df_export['percent hetatm'],df_export['hetatm 1'],df_export['hetatm 2'],df_export['hetatm 3'] = zip(*added_data)
     #print(df_export)
     os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
     df_export.to_csv('phi_xray_chains_sorted_curve_bfac_is_pos_cluster_conf_hetatmstats_mse_edo_unx__25clust_trig_SEQS_rev4.csv')


def split_hetstats_subcluster(length,is_pos,cluster,subcluster):
     os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
     phipsidf = pd.read_csv("phi_xray_chains_sorted_curve_bfac_is_pos_cluster_subcluster_conf_hetatmmulti_mse_edo_unx_25clust_TRIG_rev4.csv",usecols = ['line_number','name_chain_res','length','bfac_within_stddev','curvature','is_pos','cluster','subcluster','conf_seq','hetatm','hetatm_simple'])
     dfL = deepcopy(phipsidf[(phipsidf['length']==length) & (phipsidf['is_pos'] ==is_pos) & (phipsidf['cluster']==cluster) & (phipsidf['subcluster']==subcluster)])
     df_blanks= deepcopy(dfL[(dfL['hetatm_simple']=='[]')])
     #print(df_blanks)
     percent_het = float(len(dfL)-len(df_blanks))/len(dfL)
     #print(percent_het)
     hetlist_full = []
     hetlist = []
     #hetlist.append(ast.literal_eval(dfL['hetatm simple']))
     for index,row in dfL.iterrows():
         hets = ast.literal_eval(row['hetatm_simple'])
         for het in hets:
            hetlist_full.append(het)
     hetlist_full.sort()
     hetlist_full = list(hetlist_full for hetlist_full,_ in it.groupby(hetlist_full)) #remove dup hetatms to fix duplicate issue in 21_1  
    #just append het[0] to hetlist for regular use
     print(hetlist_full)
     hetlist = [item[0] for item in hetlist_full]
     hetdict = {}
     for het in hetlist:
         if het in hetdict:
             hetdict[het] +=1
         else:
            hetdict[het] = 1
     #print(hetdict)
     hetdict = [(hetdict[key], key) for key in hetdict]
     hetdict.sort()
     hetdict.reverse()
     top_hets = []
     added_data = [percent_het]
     for i in range(3):         
         added_data.append([hetdict[i][1],hetdict[i][0]])
     #print(top_het)
     return added_data


def combine_hetstats_subcluster():
     df_export = pd.read_csv("phi_xray_chains_sorted_curve_std_bfac_is_pos_cluster_subcluster_conf_hetatmstats_mse_edo_unx_representative_align_25clust_trig_SEQS_rev4.csv",usecols = ['length','is_pos', 'cluster','subcluster','# entries','avg_curvature','std curve','consensus','representative','std rmsd'])
     added_data =[]
     for i in range(25):
        if i == 3:
            for j in range(3):
               added_data.append(split_hetstats_subcluster(5,0,i,j))
        elif i in [5,21,23]:
            for j in range(2):
                added_data.append(split_hetstats_subcluster(5,0,i,j))
        else:
            added_data.append(split_hetstats_subcluster(5,0,i,0))

     df_export['percent hetatm'],df_export['hetatm 1'],df_export['hetatm 2'],df_export['hetatm 3'] = zip(*added_data)
     #print(df_export)
     os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
     df_export.to_csv('phi_xray_chains_sorted_curve_std_bfac_is_pos_cluster_subcluster_conf_hetatmstats_mse_edo_unx_representative_align_25clust_trig_SEQSrev5.csv')

combine_hetstats_subcluster()
#split_hetstats(5,0,3)
