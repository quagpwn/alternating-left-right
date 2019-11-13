import os
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from copy import deepcopy
import itertools



def rama_plot(x,y,x_label_in,y_label_in,name,num_values):
    sns.set()
    plt.plot(x,y,".")
    plt.xlim(-180,180) # Sets axis limits
    plt.ylim(-180,180) # Sets axis limits
    plt.xticks(np.arange(-180.1,180.1,30)) # Sets ticks markers for x axis
    plt.yticks(np.arange(-180.1,180.1,30)) # Sets ticks makers for y axis
    plt.xlabel(x_label_in) # Adds axis label
    plt.ylabel(y_label_in) # Adds axis label
    plt.arrow(-180,0,360,0) # Creates an arrow
    plt.arrow(0,-180,0,360) # Creates an arrow    
    title = name[:-4]

    plt.title(title+" "+ str(num_values)+"("+str(int(num_values/5))+" entries)")
    
    fig = plt.gcf() # Creates a figure
    #fig.set_size_inches(7.0, 7.0) # Changes figure size
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/graphs/phi graphs trig/25 clusters of len 5 trig")
    fig.savefig(name, dpi=300)
    plt.close()

    #plt.show()

def rama_plot_test(x,y,x_label_in,y_label_in,name,num_values):
    sns.set()
    plt.plot(np.cos(np.deg2rad(x))*180,np.sin(np.deg2rad(y))*180,".")
    plt.xlim(-180,180) # Sets axis limits
    plt.ylim(-180,180) # Sets axis limits
    plt.xticks(np.arange(-180.1,180.1,30)) # Sets ticks markers for x axis
    plt.yticks(np.arange(-180.1,180.1,30)) # Sets ticks makers for y axis
    plt.xlabel(x_label_in) # Adds axis label
    plt.ylabel(y_label_in) # Adds axis label
    plt.arrow(-180,0,360,0) # Creates an arrow
    plt.arrow(0,-180,0,360) # Creates an arrow    
    title = name[:-4]

    plt.title(title+" "+ str(num_values)+"("+str(int(num_values/5))+" entries)")
    
    fig = plt.gcf() # Creates a figure
    #fig.set_size_inches(7.0, 7.0) # Changes figure size
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/graphs/phi graphs/20 clusters of len 5")
    fig.savefig(name, dpi=300)
    plt.close()

    #plt.show()


def get_phipsi_list(line_number,namechainres,length):
    name = (namechainres)[0:4]
    chain = (namechainres)[10:11]
    res= (namechainres)[15:]
    line_number = line_number-1
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/phipsi_xray_chains")
    filename = "pdb"+name+"_chain"+chain+"_biopython.tsv"
    angledf = pd.read_csv(filename, sep='\t',usecols=[0,1,2], names=['name_chain_res','phi','psi'])
    phi_list = []
    psi_list = []
    for index, row in angledf.iloc[line_number:line_number+length].iterrows():
        phi_list.append(row['phi'])
        psi_list.append(row['psi'])
    #angle_list = np.array(angle_list)
    return phi_list,psi_list

def prep_points(length,is_pos,cluster,name):
     os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
     phipsidf = pd.read_csv("phi_xray_chains_sorted_curve_bfac_is_pos_cluster_conf_hetatm_25clust_TRIG_rev4.csv",usecols = ['line_number','name_chain_res','length','bfac_within_stddev','curvature','is_pos','cluster','hetatm'])
     dfL = deepcopy(phipsidf[(phipsidf['length']==length) & (phipsidf['is_pos'] ==is_pos) & (phipsidf['cluster']==cluster)])
     df_data = pd.DataFrame()
     #df_data['phi','psi'] = dfL.apply(lambda row: get_phipsi_list(row['line_number'], row['name_chain_res'], row['length']),axis = 1)
     phi = []
     psi = []
     for index,row in dfL.iterrows():
         phipsilist = get_phipsi_list(row['line_number'], row['name_chain_res'], row['length'])
         phi.extend(phipsilist[0])
         psi.extend(phipsilist[1])
     df_data['phi'] = phi
     df_data['psi'] = psi
     num_values = len(df_data.index)
     #print(phi,psi)
     #print(df_data)
     rama_plot(phi,psi,'phi','psi',name,num_values)
     #rama_plot_test(phi,psi,'phi','psi',name,num_values)


def runner():
    for i in range(25):
        prep_points(5,0,i,'rama_phi_50'+str(i)+'_25clusters_trig.png')

#runner()


def rama_plot_subcluster(xs,ys,x_label_in,y_label_in,name,num_values):
    colors = itertools.cycle(["r", "g", "b",'c','m'])
    sns.set()
    
    plt.xlim(-180,180) # Sets axis limits
    plt.ylim(-180,180) # Sets axis limits
    plt.xticks(np.arange(-180.1,180.1,30)) # Sets ticks markers for x axis
    plt.yticks(np.arange(-180.1,180.1,30)) # Sets ticks makers for y axis
    plt.xlabel(x_label_in) # Adds axis label
    plt.ylabel(y_label_in) # Adds axis label
    plt.arrow(-180,0,360,0) # Creates an arrow
    plt.arrow(0,-180,0,360) # Creates an arrow    
    title = name[:-4]

    plt.title(title+" "+ str(num_values)+"("+str(int(num_values/5))+" entries)")
    for x,y in zip(xs, ys):
        plt.scatter(x,y,color = next(colors),s=4)

    fig = plt.gcf() # Creates a figure
    #fig.set_size_inches(7.0, 7.0) # Changes figure size
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/graphs/phi graphs trig/subcluster_graphs")
    fig.savefig(name, dpi=300)
    plt.close()

    #plt.show()


def get_phipsi_list_subcluster(line_number,namechainres,length):
    name = (namechainres)[0:4]
    chain = (namechainres)[10:11]
    res= (namechainres)[15:]
    line_number = line_number-1
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/phipsi_xray_chains")
    filename = "pdb"+name+"_chain"+chain+"_biopython.tsv"
    angledf = pd.read_csv(filename, sep='\t',usecols=[0,1,2], names=['name_chain_res','phi','psi'])
    phi_list = []
    psi_list = []
    for index, row in angledf.iloc[line_number:line_number+length].iterrows():
        phi_list.append(row['phi'])
        psi_list.append(row['psi'])
    #angle_list = np.array(angle_list)
    return phi_list,psi_list


def prep_points_subcluster(length, is_pos, cluster, subcluster):
     os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
     phipsidf = pd.read_csv("phi_xray_chains_sorted_curve_bfac_is_pos_cluster_subcluster_conf_hetatmmulti_mse_edo_unx_25clust_TRIG_rev4.csv",usecols = ['line_number','name_chain_res','length','bfac_within_stddev','curvature','is_pos','cluster','subcluster','conf_seq','hetatm','hetatm_simple'])
     dfL = deepcopy(phipsidf[(phipsidf['length']==length) & (phipsidf['is_pos'] ==is_pos) & (phipsidf['cluster']==cluster) & (phipsidf['subcluster'] == subcluster)])
     df_data = pd.DataFrame()
     #df_data['phi','psi'] = dfL.apply(lambda row: get_phipsi_list(row['line_number'], row['name_chain_res'], row['length']),axis = 1)
     phi = []
     psi = []
     for index,row in dfL.iterrows():
         phipsilist = get_phipsi_list_subcluster(row['line_number'], row['name_chain_res'], row['length'])
         phi.extend(phipsilist[0])
         psi.extend(phipsilist[1])
     df_data['phi'] = phi
     df_data['psi'] = psi
     num_values = len(df_data.index)
     #print(phi,psi)
     #print(df_data)
     name = str(length)+str(is_pos)+str(cluster)+'_'+str(subcluster)+"_subcluster_graph.png"
     rama_plot_subcluster(phi,psi,'phi','psi',name, num_values)
     #rama_plot_test(phi,psi,'phi','psi',name,num_values)


def run_graph_subcluster():
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    df_export = pd.read_csv("phi_xray_chains_sorted_curve_std_bfac_is_pos_cluster_subcluster_conf_hetatmstats_mse_edo_unx_representative2.5_align_25clust_trig_SEQS_rev4.csv",usecols = ['length','is_pos', 'cluster','subcluster','# entries','avg_curvature','consensus','percent hetatm','hetatm 1','hetatm 2','hetatm 3'])
    df_export.apply(lambda row: prep_points_subcluster(row['length'],row['is_pos'],row['cluster'],row['subcluster']),axis = 1)
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")

run_graph_subcluster()

#prep_points(6,0,1,'rama_phipsi_601_overlay_separate.png')
#prep_points(6,0,0,'rama_phipsi_600_overlay_separate.png')
#plt.close()
#prep_points(5,1,1,'rama_phipsi_511.png')
#prep_points(5,1,0,'rama_phipsi_510.png')
#plt.close()
#prep_points(5,0,3,'rama_phipsi_503_overlay_separate.png')
#prep_points(5,0,2,'rama_phipsi_502_overlay_separate.png')
#prep_points(5,0,1,'rama_phipsi_501_overlay_separate.png')
#prep_points(5,0,0,'rama_phipsi_500_overlay_separate.png')
#plt.close()
#prep_points(4,1,2,'rama_phipsi_412_overlay_separate.png')
#prep_points(4,1,1,'rama_phipsi_411_overlay_separate.png')
#prep_points(4,1,0,'rama_phipsi_410_overlay_separate.png')
#plt.close()
#prep_points(4,0,2,'rama_phipsi_402_overlay_separate.png')
#prep_points(4,0,1,'rama_phipsi_401_overlay_separate.png')
#prep_points(4,0,0,'rama_phipsi_400_overlay_separate.png')
#plt.close()

#prep_points(11,0,1,'rama_phi_1101.png')
#prep_points(11,0,0,'rama_phi_1100.png')
#plt.close()
#prep_points(9,0,1,'rama_phi_901.png')
#prep_points(9,0,0,'rama_phi_900.png')
#plt.close()
#prep_points(8,1,1,'rama_phi_811.png')
#prep_points(8,1,0,'rama_phi_810.png')
#plt.close()
#prep_points(8,0,1,'rama_phi_801.png')
#prep_points(8,0,0,'rama_phi_800.png')
#plt.close()
#prep_points(7,0,2,'rama_phi_702.png')
#prep_points(7,0,1,'rama_phi_701.png')
#prep_points(7,0,0,'rama_phi_700.png')
#plt.close()
#prep_points(6,1,5,'rama_phi_615.png')
#prep_points(6,1,4,'rama_phi_614.png')
#prep_points(6,1,3,'rama_phi_613.png')
#prep_points(6,1,2,'rama_phi_612.png')
#prep_points(6,1,1,'rama_phi_611.png')
#prep_points(6,1,0,'rama_phi_610.png')
#plt.close()
#prep_points(6,0,3,'rama_phi_603.png')
#prep_points(6,0,2,'rama_phi_602.png')
#prep_points(6,0,1,'rama_phi_601.png')
#prep_points(6,0,0,'rama_phi_600.png')
#plt.close()
#prep_points(5,1,1,'rama_phi_511.png')
#prep_points(5,1,0,'rama_phi_510.png')
#plt.close()
#prep_points(4,1,4,'rama_phi_414.png')
#prep_points(4,1,3,'rama_phi_413.png')
#prep_points(4,1,2,'rama_phi_412.png')
#prep_points(4,1,1,'rama_phi_411.png')
#prep_points(4,1,0,'rama_phi_410.png')
#plt.close()
#prep_points(4,0,4,'rama_phi_404.png')
#prep_points(4,0,3,'rama_phi_403.png')
#prep_points(4,0,2,'rama_phi_402.png')
#prep_points(4,0,1,'rama_phi_401.png')
#prep_points(4,0,0,'rama_phi_400.png')
#plt.close()


#prep_points(5,0,19,'rama_phi_5019_20clusters_trig.png')
#prep_points(5,0,18,'rama_phi_5018_20clusters_trig.png')
#prep_points(5,0,17,'rama_phi_5017_20clusters_trig.png')
#prep_points(5,0,16,'rama_phi_5016_20clusters_trig.png')
#prep_points(5,0,15,'rama_phi_5019_20clusters_trig.png')
#prep_points(5,0,14,'rama_phi_5014_20clusters_trig.png')
#prep_points(5,0,13,'rama_phi_5013_20clusters_trig.png')
#prep_points(5,0,12,'rama_phi_5012_20clusters_trig.png')
#prep_points(5,0,11,'rama_phi_5011_20clusters_trig.png')
#prep_points(5,0,10,'rama_phi_5010_20clusters_trig.png')
#prep_points(5,0,9,'rama_phi_509_20clusters_trig.png')
#prep_points(5,0,8,'rama_phi_508_20clusters_trig.png')
#prep_points(5,0,7,'rama_phi_507_20clusters_trig.png')
#prep_points(5,0,6,'rama_phi_506_20clusters_trig.png')
#prep_points(5,0,5,'rama_phi_505_20clusters_trig.png')
#prep_points(5,0,4,'rama_phi_504_20clusters_trig.png')
#prep_points(5,0,3,'rama_phi_503_20clusters_trig.png')
#prep_points(5,0,2,'rama_phi_502_20clusters_trig.png')
#prep_points(5,0,1,'rama_phi_501_20clusters_trig.png')
#prep_points(5,0,0,'rama_phi_500_20clusters_trig.png')






