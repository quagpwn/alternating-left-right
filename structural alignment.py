import os
import time
import numpy as np
#import math
import sys
from Bio.PDB import *
import pandas as pd
from copy import deepcopy
import itertools
import __main__
os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/repos")
from align_allfiles_to_allfiles import *
import collections as ct
import pprint
#import tmalign

__main__.pymol_argv = ['pymol','-c'] # Pymol: quiet and no GUI
import pymol
pymol.finish_launching()
print('yuh')

#pymol.cmd.load(mobileStructurePath, mobileStructureName)
        
def multi_align(length,is_pos,cluster,subcluster):
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    phipsidf = pd.read_csv("phi_xray_chains_sorted_curve_bfac_is_pos_cluster_subcluster_conf_hetatmmulti_mse_edo_unx_25clust_TRIG_rev4.csv",usecols = ['line_number','name_chain_res','length','bfac_within_stddev','curvature','is_pos','cluster','subcluster','conf_seq','hetatm_simple','hetatm'])
    dfL = deepcopy(phipsidf[(phipsidf['length']==length) & (phipsidf['is_pos'] ==is_pos) & (phipsidf['cluster']==cluster) & (phipsidf['subcluster']==subcluster)])
    fldr = "C:/Users/Shaheer Rizwan/Documents/ramachandran/cropped_pdb_subcluster_rev2/"+str(length)+str(is_pos)+str(cluster)+"_"+str(subcluster)
    os.chdir(fldr)
    #print(os.getcwd())
    df_name = pd.Series.to_list(dfL['name_chain_res'].str.slice(start=0,stop=4)+'_chain'+dfL['name_chain_res'].str.slice(start=10,stop=11)+'_'+dfL['name_chain_res'].str.slice(start=12)+'_cropped.pdb')
    #print(df_name)
    dd = ct.defaultdict(list)

    for ref, comp in itertools.combinations(dfL['name_chain_res'],2):
        pymol.cmd.reinitialize()
        ref_name = (ref)[0:4]
        ref_chain = (ref)[10:11]
        ref_res = (ref)[12:15]
        ref_resid= (ref)[15:]

        comp_name = comp[0:4]
        comp_chain = (comp)[10:11]
        comp_res = (comp)[12:15]
        comp_resid= (comp)[15:]

        ref_name = ref_name + '_chain'+str(ref_chain)+'_'+str(ref_res)+str(ref_resid)
        comp_name = comp_name + '_chain'+str(comp_chain)+'_'+str(comp_res)+str(comp_resid)

        #print(ref_name,comp_name)
        pymol.cmd.load(ref_name+'_cropped.pdb',ref_name)
        pymol.cmd.load(comp_name+'_cropped.pdb',comp_name)
        #rmsd_raw=pymol.cmd.super(ref_name+" " +" c. a & n. ca", comp_name+" " +" c. a & n. ca" , quiet=0, object = ref_name+'_'+comp_name) #v1
        rmsd_raw=pymol.cmd.super(ref_name, comp_name, quiet=0, cycles = 0, object = ref_name+'_'+comp_name) #v2 in progress - trying to do no outlier rejection
        #pymol.cmd.save(filename = ref_name+'_'+comp_name+'.pse',)
        #time.sleep(1)
        dd[ref_name].append(rmsd_raw[0])
        dd[comp_name].append(rmsd_raw[0])

    dd_mean = ct.defaultdict(list)
    for ref in dd:
        #print(dd[ref])
        dd_mean[ref] = np.mean(dd[ref])
    #pprint.pprint(dd)
    print(dd_mean)
    representative = min(dd,key=dd.get)
    std = np.std(list(dd_mean.values()))
    print(length,is_pos,cluster,subcluster,representative, std)
    return representative,std


#    pymol.cmd.reinitialize()
#    for fl in df_name:
##        print('loading '+fl)
#        pymol.cmd.load(fl,fl[:-4])
#    #print(df_name[0][0:4])
#    df_selects = []
#    rmsds = []
#    #all_aligned = pymol.cmd.extra_fit(reference = representative+" " +"& backbone", method = 'super',quiet = 0,object = 'alignment_object')
#    for comp in dfL['name_chain_res']:
#        comp_name = comp[0:4]
        
#        rmsd_raw = pymol.cmd.super(comp_name, '1kqf' , quiet=0)
#        print(comp_name+ " "+str(rmsd_raw[0]))
#        rmsds.append(rmsd_raw[0])
        
#    rmsd_cutoff = np.percentile(rmsds,25)
#    print(rmsd_cutoff)

#     #start_time = time.time()
#    pymol.cmd.reinitialize()
#    for fl in df_name:
##        print('loading '+fl)
#        pymol.cmd.load(fl,fl[:-4])
#    #print(df_name[0][0:4])

#    for  comp in dfL['name_chain_res']:  
#        comp_name = comp[0:4]
        
#        rmsd_raw = pymol.cmd.super(comp_name, '1kqf' , quiet=0)
#        if rmsd_raw[0] < rmsd_cutoff:
#            df_selects.append(comp)
#        #print(comp_name+ " "+str(rmsd_raw[0]))
#    print(df_selects)



#    pymol.cmd.reinitialize()
#    for fl in df_selects:
##        print('loading '+fl)
#        pymol.cmd.load(fl[0:4]+'_cropped.pdb',fl[0:4])
#    #print(df_name[0][0:4])

#    for comp in df_selects:  
#        comp_name = comp[0:4]        
#        pymol.cmd.super(comp_name+' & backbone', '1kqf & backbone' , quiet=0 )
        
#    print(len(df_name))
#    print(len(df_selects))
#    pymol.cmd.save(filename = str(length)+str(is_pos)+str(cluster)+'_all.pse')
    #print(all_aligned)
    #print(representative)
    ##print(rmsd_list)
    #end_time = time.time()
    #print (str(length)+str(is_pos)+str(cluster)+' Run time is:'+str(end_time - start_time))
    #return np.mean(rmsd_list)



def biopython_multi_align(length,is_pos,cluster): #biopython alignment doesnt work
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    phipsidf = pd.read_csv("phi_xray_chains_sorted_curve_bfac_is_pos_cluster_conf_hetatmsimple_mse_25clust_TRIG_rev4.csv",usecols = ['line_number','name_chain_res','length','bfac_within_stddev','curvature','is_pos','cluster','conf_seq','hetatm simple','hetatm'])
    dfL = deepcopy(phipsidf[(phipsidf['length']==length) & (phipsidf['is_pos'] ==is_pos) & (phipsidf['cluster']==cluster)])
    fldr = "C:/Users/Shaheer Rizwan/Documents/ramachandran/cropped_pdb/"+str(length)+str(is_pos)+str(cluster)
    os.chdir(fldr)
    print(os.getcwd())
    df_name = pd.Series.to_list(dfL['name_chain_res'].str.slice(start=0,stop=4)+'_cropped.pdb')
    print(df_name)
    for ref,comp in itertools.combinations(dfL['name_chain_res'],2):
        ref_name = (ref)[0:4]
        comp_name = comp[0:4]
        ref_structure = PDBParser().get_structure(ref_name, ref_name+'_cropped.pdb')
        comp_structure = PDBParser().get_structure(comp_name, comp_name+'_cropped.pdb')
        ref_model = ref_structure[0]
        comp_model = comp_structure[0]
        ref_atoms = []
        comp_atoms = []

        for chain in ref_model:
            for res in chain:
                for atoms in res:
                    ref_atoms.append(atoms)

        for chain in comp_model:
            for res in chain:
                for atoms in res:
                    comp_atoms.append(atoms)

        super_imposer = Superimposer()
        super_imposer.set_atoms(ref_atoms, comp_atoms)
        super_imposer.apply(comp_atoms)
        print(super_imposer.rms)

        io = Bio.PDB.PDBIO()
        io.set_structure(sample_structure) 
        io.save(ref_name+'_'+comp_name+".pdb")


def align_clusters():
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    df_export = pd.read_csv("phi_xray_chains_sorted_curve_bfac_is_pos_cluster_conf_hetatmstats_mse_edo_unx__25clust_trig_SEQS_rev4.csv",usecols = ['length','is_pos', 'cluster','# entries','avg_curvature','consensus','percent hetatm','hetatm 1','hetatm 2','hetatm 3'])
    df_export['avg_rmsd_all'] = df_export.apply(lambda row: multi_align(row['length'],row['is_pos'],row['cluster']),axis = 1)
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    df_export.to_csv('phi_xray_chains_sorted_curve_bfac_is_pos_cluster_conf_hetatmstats_mse_edo_unx__align_25clust_trig_SEQS_rev4.csv')
     

def write_representative():
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    df_export = pd.read_csv("phi_xray_chains_sorted_curve_std_bfac_is_pos_cluster_subcluster_conf_hetatmstats_mse_edo_unx_representative2_align_25clust_trig_SEQS_rev4.csv",usecols = ['length','is_pos', 'cluster','subcluster','# entries','avg_curvature','std curve','consensus','percent hetatm','hetatm 1','hetatm 2','hetatm 3'])
    df_export[['representative','std rmsd']] = df_export.apply(lambda row: multi_align(row['length'],row['is_pos'],row['cluster'],row['subcluster']),axis = 1,result_type = 'expand')
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    df_export.to_csv('phi_xray_chains_sorted_curve_std_bfac_is_pos_cluster_subcluster_conf_hetatmstats_mse_edo_unx_representative2.5_align_25clust_trig_SEQS_rev4.csv')

start = time.time()
#print(multi_align(5,0,23,0))

write_representative()
end = time.time()
print ('Run time is:'+str(end - start))
pymol.cmd.quit()
