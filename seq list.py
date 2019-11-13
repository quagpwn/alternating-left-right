from Bio.PDB import *
import os
import pandas as pd
import numpy as np
import math
import os
#import weblogolib as w
import io
from Bio import motifs
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from copy import deepcopy
#import __main__
#__main__.pymol_argv = ['pymol','-c'] # Pymol: quiet and no GUI
#import pymol
#pymol.finish_launching()
#print('yuh')


def get_seqlist(length,is_pos,cluster,subcluster):
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    phipsidf = pd.read_csv("phi_xray_chains_sorted_curve_bfac_is_pos_cluster_subcluster_conf_hetatmmulti_mse_edo_unx_25clust_TRIG_rev4.csv",usecols = ['line_number','name_chain_res','length','bfac_within_stddev','curvature','is_pos','cluster','subcluster','conf_seq','hetatm_simple','hetatm'])
    dfL = deepcopy(phipsidf[(phipsidf['length']==length) & (phipsidf['is_pos'] ==is_pos) & (phipsidf['cluster']==cluster) & (phipsidf['subcluster']==subcluster)])
    fldr = "C:/Users/Shaheer Rizwan/Documents/ramachandran/cropped_pdb_subcluster_rev2/"+str(length)+str(is_pos)+str(cluster)+"_"+str(subcluster)
    os.chdir(fldr)
    #print(os.getcwd())
    df_name = pd.Series.to_list(dfL['name_chain_res'].str.slice(start=0,stop=4)+'_cropped.pdb')
    seq_dict = {}
    for ref in dfL['name_chain_res']:
        #pymol.cmd.reinitialize()
        ref_name = (ref)[0:4]
        ref_chain = (ref)[10:11]
        ref_res = (ref)[12:15]
        ref_resid= (ref)[15:]
        ref_name = ref_name + '_chain'+str(ref_chain)+'_'+str(ref_res)+str(ref_resid)
        ref_name_pdb = ref_name+'_cropped.pdb'
        #pymol.cmd.load(ref_name+'_cropped.pdb',ref_name)
        #seq = cmd.get_fastastr('all') get seq if using pymol
        structure = PDBParser().get_structure(ref_name,ref_name_pdb)
        ppb = CaPPBuilder()
        for pp in ppb.build_peptides(structure):
            seq = pp.get_sequence()
            print(seq, str(seq))
            seq_dict[ref_name] = str(seq)
    seq_list = list(seq_dict.values())
    df_seq = pd.DataFrame(seq_list)
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/sequence_lists_rev2")
    df_seq.to_csv(str(length)+str(is_pos)+str(cluster)+"_"+str(subcluster)+"_sequences.csv")
    return seq_list

def compile_seqtable():
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    df_export = pd.read_csv("phi_xray_chains_sorted_curve_std_bfac_is_pos_cluster_subcluster_conf_hetatmstats_mse_edo_unx_representative2_align_25clust_trig_SEQS_rev4.csv",usecols = ['length','is_pos', 'cluster','subcluster','# entries','avg_curvature','consensus','percent hetatm','hetatm 1','hetatm 2','hetatm 3'])
    df_seqs = pd.DataFrame()
    df_seqs['seq list'] = df_export.apply(lambda row: get_seqlist(row['length'],row['is_pos'],row['cluster'],row['subcluster']),axis = 1)
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    #df_export.to_csv('phi_xray_chains_sorted_curve_bfac_is_pos_cluster_subcluster_conf_hetatmstats_mse_edo_unx_representative_align_25clust_trig_SEQS_seqs_rev4.csv')


compile_seqtable()