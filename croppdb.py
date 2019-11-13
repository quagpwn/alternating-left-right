import os
import pandas as pd
import numpy as np
import pymol
from copy import deepcopy
import ast
from Bio import PDB
from Bio.PDB import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import *
import time




def get_cropped_pdb(namechainres,length,is_pos,cluster,hetlist):
    name = (namechainres)[0:4]
    chain = (namechainres)[10:11]
    print(namechainres)
    resid= (namechainres)[15:]
    #print(type(resid))
    resid = int(resid)
    
    length = int(length)    
    hetlist_full = ast.literal_eval(hetlist)
    hetlist = []
    for subl in hetlist_full:
        hetlist.append(subl[1])
        
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/ids_xray")
    io = PDB.PDBIO()
    structure = PDBParser().get_structure(name,'pdb'+name+".ent")
    class HeterochiralSelect(Select):
        def accept_residue(self, residue):
            #print(residue.id)
            if ((residue.id[1] in range(resid,length+resid)) or residue.id[1] in hetlist):
                #print(residue['CA'].get_serial_number())
                return 1
            else:
                return 0

    io.set_structure(structure[0][chain])
    if not os.path.exists("C:/Users/Shaheer Rizwan/Documents/ramachandran/cropped_pdb/"+str(length)+str(is_pos)+str(cluster)):
        os.makedirs("C:/Users/Shaheer Rizwan/Documents/ramachandran/cropped_pdb/"+str(length)+str(is_pos)+str(cluster))
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/cropped_pdb/"+str(length)+str(is_pos)+str(cluster))
    io.save(name + '_cropped.pdb', HeterochiralSelect())

def get_cropped_pdb_no_het(namechainres,length,is_pos,cluster,hetlist):
    name = (namechainres)[0:4]
    chain = (namechainres)[10:11]
    print(namechainres)
    resid= (namechainres)[15:]
    #print(type(resid))
    resid = int(resid)
    
    length = int(length)    
    #hetlist_full = ast.literal_eval(hetlist)
    #hetlist = []
    #for subl in hetlist_full:
    #    hetlist.append(subl[1])
        
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/ids_xray")
    io = PDB.PDBIO()
    structure = PDBParser().get_structure(name,'pdb'+name+".ent")
    class HeterochiralSelect(Select):
        def accept_residue(self, residue):
            #print(residue.id)
            if (residue.id[1] in range(resid,length+resid)):
                #print(residue['CA'].get_serial_number())
                return 1
            else:
                return 0

    io.set_structure(structure[0][chain])
    if not os.path.exists("C:/Users/Shaheer Rizwan/Documents/ramachandran/cropped_pdb_no_het/"+str(length)+str(is_pos)+str(cluster)):
        os.makedirs("C:/Users/Shaheer Rizwan/Documents/ramachandran/cropped_pdb_no_het/"+str(length)+str(is_pos)+str(cluster))
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/cropped_pdb_no_het/"+str(length)+str(is_pos)+str(cluster))
    io.save(name + '_cropped_no_het.pdb', HeterochiralSelect())

def get_cropped_pdb_subcluster(namechainres,length,is_pos,cluster,subcluster,hetlist):
    name = (namechainres)[0:4]
    chain = (namechainres)[10:11]
    #print(namechainres)
    res = (namechainres)[12:15]
    resid= (namechainres)[15:]
    #print(type(resid))
    resid = int(resid)
    length = int(length)
    print(name + '_chain'+str(chain)+'_'+str(res)+str(resid)+'_cropped.pdb')

    hetlist_full = ast.literal_eval(hetlist)
    hetlist = []
    for subl in hetlist_full:
        hetlist.append(subl[1])
        
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/ids_xray")
    io = PDB.PDBIO()
    structure = PDBParser().get_structure(name,'pdb'+name+".ent")
    class HeterochiralSelect(Select):
        def accept_residue(self, residue):
            #print(residue.id)
            if ((residue.id[1] in range(resid,length+resid)) or residue.id[1] in hetlist):
                #print(residue['CA'].get_serial_number())
                return 1
            else:
                return 0

    io.set_structure(structure[0][chain])
    if not os.path.exists("C:/Users/Shaheer Rizwan/Documents/ramachandran/cropped_pdb_subcluster_rev2/"+str(length)+str(is_pos)+str(cluster)+"_"+str(subcluster)):
        os.makedirs("C:/Users/Shaheer Rizwan/Documents/ramachandran/cropped_pdb_subcluster_rev2/"+str(length)+str(is_pos)+str(cluster)+"_"+str(subcluster))
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/cropped_pdb_subcluster_rev2/"+str(length)+str(is_pos)+str(cluster)+"_"+str(subcluster))
    io.save(name + '_chain'+str(chain)+'_'+str(res)+str(resid)+'_cropped.pdb', HeterochiralSelect())



def hetcalc(namechainres,length):
    name = (namechainres)[0:4]
    chain = (namechainres)[10:11]
    print(namechainres)
    resid= (namechainres)[15:]
    resid = int(resid)
    #print(resid)
    length = int(length)    
    hetlist =[]
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/ids_xray")
    structure = PDBParser().get_structure(name,'pdb'+name+".ent")
    for residue in range(resid, resid+length):
        try:
            nitrogen = structure[0][chain][residue]['N']
            #print(nitrogen.get_parent().get_resname(),residue)
            for hetres in structure[0][chain]:
                if hetres.id[0] != ' ':
                    for hetatm in hetres:
                        #print(hetatm.get_name())
                        dist = (nitrogen - hetatm)
                        if dist < 5.0 and hetatm.get_parent().get_resname() != "HOH" and hetatm.get_parent().get_resname() != "MSE":
                            hetlist.append([nitrogen.get_parent().get_resname()+str(residue),hetatm.get_name(),hetatm.get_parent().get_resname(), str(hetatm.get_serial_number())])
                            #print(dist,hetatm.get_name(),hetatm.get_parent().get_resname(), hetatm.get_serial_number(),"\n")
        except KeyError:
            continue
    return hetlist


def crop_pdb_list(length):
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    phipsidf = pd.read_csv("phi_xray_chains_sorted_curve_bfac_is_pos_cluster_subcluster_conf_hetatmmulti_mse_edo_unx_25clust_TRIG_rev4.csv",usecols = ['line_number','name_chain_res','length','bfac_within_stddev','curvature','is_pos','cluster','subcluster','conf_seq','hetatm','hetatm_simple'])
    dfL = deepcopy(phipsidf[(phipsidf['length']==length)])
    dfL.apply(lambda row: get_cropped_pdb_subcluster(row['name_chain_res'], row['length'],row['is_pos'],row['cluster'],row['subcluster'],row['hetatm_simple']),axis=1)


start = time.time()
#get_cropped_pdb_no_het('4l8a:ChainA:ARG88',5,0,3,"[['COA', 204], ['SO4', 206]]")
crop_pdb_list(5)
end = time.time()
print('Run time is:'+str(end - start))