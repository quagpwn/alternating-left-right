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

#class RefSeqColor(w.ColorRule): #unused
#    """
#    Color the given reference sequence in its own color, so you can easily see 
#    which positions match that sequence and which don't.
#    """

#    def __init__(self, ref_seq, color, description=None):
#        self.ref_seq = ref_seq
#        self.color = w.Color.from_string(color)
#        self.description = description

#    def symbol_color(self, seq_index, symbol, rank):
#        if symbol == self.ref_seq[seq_index]:
#            return self.color

#baserules = [
#            w.SymbolColor("PBGSTYC", "green", "polar"),
#            w.SymbolColor("NQ", "purple", "neutral"),
#            w.SymbolColor("KRH", "blue", "basic"),
#            w.SymbolColor("DEL", "red", "acidic"),
#            w.SymbolColor("AWFIMV", "black", "hydrophobic")
            
#        ]

#protein_alphabet = w.Alphabet('ACDEFGHIKLMNOPQRSTUVWYBJZX*-adefghiklmnopqrstuvwybjzx', [])

def plotseqlogo(refseq, mseqs, name): #unused
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    fasta = "> \n" + "\n> \n".join(mseqs)
    seqs = w.read_seq_data(io.StringIO(fasta), alphabet=protein_alphabet)

    colorscheme = w.ColorScheme([RefSeqColor(refseq, "orange", "refseq")] + baserules,
                                alphabet = protein_alphabet)

    data = w.LogoData.from_seqs(seqs)
    options = w.LogoOptions()
    # options.logo_title = name
    options.show_fineprint = False
    options.yaxis_label = ""
    options.color_scheme = colorscheme
    options.title = name
    options.resolution=400
    mformat = w.LogoFormat(data, options)
    mimage = w.pdf_formatter(data,mformat)
    
    fname = "%s.pdf" % name
    f = open(fname, "w")
    f.write(mimage)
    f.close()


def get_phipsi_bin(line_number,namechainres,length,is_pos,cluster):
    name = (namechainres)[0:4]
    chain = (namechainres)[10:11]
    res= (namechainres)[15:]
    line_number = line_number-1
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/phipsi_xray_chains")
    filename = "pdb"+name+"_chain"+chain+"_biopython.tsv"
    angledf = pd.read_csv(filename, sep='\t',usecols=[1,2], names=['phi','psi'])
    angledf = angledf.iloc[line_number:line_number+length]
    #angledf.drop(['name_chain_res'],axis=1)
    angledf['bin'] = angledf.apply(lambda row: bin_angles(row['phi'],row['psi']), axis = 1)
    binseq = ''
    for index,row in angledf.iterrows():
        binseq+=row['bin']

    return binseq

def bin_angles(phi,psi):
    binseq = ''
    if (((phi < -30.0) and (phi > -150.0)) and ((psi < 45.0) and (psi > -90.0))): 
        binseq+='aR'
    elif (((phi < 150.0) & (phi > 30.0)) & ((psi < 90.0) & (psi > -45.0))):
        binseq+='aL'
    elif (((phi <-90.0) & (phi > -150.0)) & (((psi > 90.0) & (psi < 180.0)) | ((psi > -180.0) & (psi < -120.0)))):
        binseq+='bR'
    elif (((phi > 90.0) & (phi < 150.0)) & (((psi > -180.0) & (psi < -90.0)) | ((psi > 120.0) & (psi < -180.0)))):
        binseq+='bL'
    elif (((phi > -90.0) & (phi < -30.0)) & (((psi > 90.0) & (psi < 180.0)) | ((psi > -180.0) & (psi < -120.0)))):
        binseq+='pR'
    elif(((phi >30.0) & (phi < 90.0)) & (((psi > -180.0) & (psi < -90.0)) | ((psi > 120.0) & (psi < -180.0)))):
        binseq+='pL'
    else:
        binseq+='NA'
    return binseq



def run_bin():
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    phipsidf = pd.read_csv("phi_xray_chains_sorted_curve_bfac_is_pos_cluster_conf_hetatmsimple_mse_25clust_TRIG_rev4.csv",usecols = ['line_number','name_chain_res','length','bfac_within_stddev','curvature','is_pos','cluster','hetatm'])
    #phipsidf = phipsidf.sort_values(by ='length', ascending = False)
    phipsidf['conf_seq'] = phipsidf.apply(lambda row: get_phipsi_bin(row['line_number'],row['name_chain_res'],row['length'],row['is_pos'],row['cluster']), axis = 1)
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    phipsidf.to_csv('phi_xray_chains_sorted_curve_bfac_is_pos_cluster_conf_hetatmsimple_mse_25clust_TRIG_rev4.csv')
    return phipsidf

def run_consensus_cluster(length,is_pos,cluster):
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    phipsidf = pd.read_csv("phi_xray_chains_sorted_curve_bfac_is_pos_cluster_conf_hetatmsimple_mse_25clust_TRIG_rev4.csv",usecols = ['line_number','name_chain_res','length','bfac_within_stddev','curvature','is_pos','cluster','hetatm','conf_seq'])
    #phipsidf = phipsidf.sort_values(by ='length', ascending = False)
    dfL = deepcopy(phipsidf[(phipsidf['length']==length) & (phipsidf['is_pos'] ==is_pos) & (phipsidf['cluster']==cluster)])
    #print(dfL)    
    df_instances =  pd.DataFrame()
    df_instances['seqs'] = dfL.apply(lambda row: Seq(row['conf_seq'].upper(), IUPAC.ExtendedIUPACProtein),axis = 1)
    df_string_instances = pd.DataFrame()
    df_string_instances = df_instances.apply(lambda row: str(row['seqs']),axis = 1)
    #print(df_instances)
    m = motifs.create(df_instances['seqs'])
    consensus_dict = {'length':length,'is_pos':is_pos,'cluster':cluster,'avg_curvature': dfL['curvature'].mean(), 'consensus':str(m.consensus), '# entries': len(dfL)}
    #return consensus_dict
    consensusdf = pd.DataFrame.from_records(consensus_dict,columns = ['length','is_pos','cluster','avg_curvature','consensus', '# entries'],index = [0]) 
    return consensusdf
    #plotseqlogo(m.consensus,df_string_instances,'bigboy')
    #print(m.consensus)
    #os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    ##phipsidf.to_csv('phi_xray_chains_sorted_curve_bfac_is_pos_cluster_conf_hetatm_50clust_TRIG_rev4.csv')

def runner():
    df_compare = pd.DataFrame(columns = ['length','is_pos','cluster','avg_curvature','consensus'])
    for i in range(25):
        df_compare = df_compare.append(run_consensus_cluster(5,0,i))
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    df_compare.to_csv('phi_xray_chains_sorted_curve_bfac_is_pos_cluster_conf_hetatmstats_mse_25clust_trig_SEQS_rev4.csv')


def run_consensus_subcluster(length,is_pos,cluster,subcluster):
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    phipsidf = pd.read_csv("phi_xray_chains_sorted_curve_bfac_is_pos_cluster_subcluster_conf_hetatmmulti_mse_edo_unx_25clust_TRIG_rev4.csv",usecols = ['line_number','name_chain_res','length','bfac_within_stddev','curvature','is_pos','cluster','subcluster','conf_seq','hetatm','hetatm_simple'])
    #phipsidf = phipsidf.sort_values(by ='length', ascending = False)
    dfL = deepcopy(phipsidf[(phipsidf['length']==length) & (phipsidf['is_pos'] ==is_pos) & (phipsidf['cluster']==cluster) & (phipsidf['subcluster']==subcluster)])
    #print(dfL)    
    df_instances =  pd.DataFrame()
    df_instances['seqs'] = dfL.apply(lambda row: Seq(row['conf_seq'].upper(), IUPAC.ExtendedIUPACProtein),axis = 1)
    df_string_instances = pd.DataFrame()
    df_string_instances = df_instances.apply(lambda row: str(row['seqs']),axis = 1)
    #print(df_instances)
    m = motifs.create(df_instances['seqs'])
    consensus_dict = {'length':length,'is_pos':is_pos,'cluster':cluster, 'subcluster':subcluster, 'avg_curvature': dfL['curvature'].mean(), 'consensus':str(m.consensus), '# entries': len(dfL)}
    #return consensus_dict
    consensusdf = pd.DataFrame.from_records(consensus_dict,columns = ['length','is_pos','cluster','subcluster','avg_curvature','consensus', '# entries'],index = [0]) 
    return consensusdf
    #plotseqlogo(m.consensus,df_string_instances,'bigboy')
    #print(m.consensus)
    #os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    ##phipsidf.to_csv('phi_xray_chains_sorted_curve_bfac_is_pos_cluster_conf_hetatm_50clust_TRIG_rev4.csv')

def runner_subcluster():
    df_compare = pd.DataFrame(columns = ['length','is_pos','cluster','subcluster','avg_curvature','consensus','# entries'])
    for i in range(25):
        if i == 3:
            for j in range(3):
                df_compare = df_compare.append(run_consensus_subcluster(5,0,i,j))
        elif i in [5,21,23]:
            for j in range(2):
                df_compare = df_compare.append(run_consensus_subcluster(5,0,i,j))
        else:
            df_compare = df_compare.append(run_consensus_subcluster(5,0,i,0))
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    df_compare.to_csv('phi_xray_chains_sorted_curve_bfac_is_pos_cluster_subcluster_conf_hetatmstats_mse_edo_unx__align_25clust_trig_SEQS_rev4.csv')




def add_std_curvature(length,is_pos,cluster,subcluster):
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    phipsidf = pd.read_csv("phi_xray_chains_sorted_curve_bfac_is_pos_cluster_subcluster_conf_hetatmmulti_mse_edo_unx_25clust_TRIG_rev4.csv",usecols = ['line_number','name_chain_res','length','bfac_within_stddev','curvature','is_pos','cluster','subcluster','conf_seq','hetatm','hetatm_simple'])
    #phipsidf = phipsidf.sort_values(by ='length', ascending = False)
    dfL = deepcopy(phipsidf[(phipsidf['length']==length) & (phipsidf['is_pos'] ==is_pos) & (phipsidf['cluster']==cluster) & (phipsidf['subcluster']==subcluster)])
    return np.std(dfL['curvature'])

def runner_std_curvature():
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    df_export = pd.read_csv("phi_xray_chains_sorted_curve_bfac_is_pos_cluster_subcluster_conf_hetatmstats_mse_edo_unx_representative_align_25clust_trig_SEQS_rev4.csv",usecols = ['length','is_pos', 'cluster','subcluster','# entries','avg_curvature','consensus','percent hetatm','hetatm 1','hetatm 2','hetatm 3','representative','std rmsd'])
    df_export['std curve'] = df_export.apply(lambda row: add_std_curvature(row['length'],row['is_pos'],row['cluster'],row['subcluster']),axis = 1)
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    df_export.to_csv('phi_xray_chains_sorted_curve_std_bfac_is_pos_cluster_subcluster_conf_hetatmstats_mse_edo_unx_representative_align_25clust_trig_SEQS_rev4.csv')
    
    
#runner_subcluster()
#run_bin()
#print(run_consensus_cluster(5,0,49
#print(get_phipsi_bin(337,"5a57:ChainA:GLY653",11))
runner_std_curvature()