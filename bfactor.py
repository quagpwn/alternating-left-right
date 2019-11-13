from Bio.PDB import *
import os
import pandas as pd
import numpy as np

def calc_bfac(namechainres,length):
    name = (namechainres)[0:4]
    chainid = (namechainres)[10:11]
    #print(namechainres)
    resid= (namechainres)[15:]
    resid = int(resid)
    #print(resid)
    length = int(length)    
    hetlist =[]
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/ids_xray")
    structure = PDBParser().get_structure(name,'pdb'+name+".ent")
    chain = structure[0][chainid]
    for residue in chain:
        for atom in residue:
            try:
                hetlist.append(atom.get_bfactor())
                #print(atom.get_bfactor())
            except KeyError:
                continue                
    b_avg = np.mean(hetlist)
    b_stddev = np.std(hetlist)
    #print('mean',b_avg)
    hetchiral_list = []
    for residue in range(resid, resid+length):
        try:
            for atom in structure[0][chainid][residue]:
                try:
                    atom_bfac = atom.get_bfactor()
                    hetchiral_list.append(atom_bfac)
                except KeyError:
                    continue
        except KeyError:
            continue
            #print(atom_bfac)
    hetchiral_avg = np.mean(hetchiral_list)
    #print("heterochiral",hetchiral_avg-b_avg,b_stddev)
    if ((hetchiral_avg-b_avg) >= b_stddev):
        return 0
    else:
        return 1


def bfac():
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    phipsidf = pd.read_csv("phi_xray_chains_sorted_curvature_is_pos_cluster_hetatm_rev4-testing5-20clusters.csv",usecols = ['line_number','name_chain_res','length','curvature','is_pos','cluster','hetatm'])
    phipsidf['bfac_within_stddev'] = phipsidf.apply(lambda row: calc_bfac(row['name_chain_res'],row['length']), axis = 1)
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    phipsidf.to_csv('phi_xray_chains_sorted_curvature_bfac_is_pos_cluster_hetatm_rev4-testing5-20clusters.csv')
    return phipsidf

#print(calc_bfac("4bz4:ChainA:TYR208",11))

print(bfac())
    






