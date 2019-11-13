from Bio.PDB import *
import os
import pandas as pd
import numpy as np
import math

def degrees(rad_angle) :
    """Converts any angle in radians to degrees.
    If the input is None, the it returns None.
    For numerical input, the output is mapped to [-180,180]
    """
    if rad_angle is None :
        return None
    angle = rad_angle * 180 / math.pi
    while angle > 180 :
        angle = angle - 360
    while angle < -180 :
        angle = angle + 360
    return angle

def calccurve(namechainres,length):
    name = (namechainres)[0:4]
    chain = (namechainres)[10:11]
    resid= (namechainres)[15:]
    resid = int(resid)
    length = int(length)    
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/ids_xray")
    structure = PDBParser().get_structure(name,'pdb'+name+".ent")
    c_vectors = []
    for residue in range(resid, resid+length,2):
        try:
            c_alpha = structure[0][chain][residue]['CA']
            c_alpha_angle = c_alpha.get_vector()
            print(c_alpha.get_parent().get_resname(),residue,c_alpha_angle)
            c_vectors.append(c_alpha_angle)

        except KeyError:
            return np.NaN  
    angle = calc_angle(c_vectors[0],c_vectors[1],c_vectors[2])
    return degrees(angle)

def calccurve_general(namechainres,length):
    if length == 5:
        return calccurve(namechainres,length)
    else:
        name = (namechainres)[0:4]
        chain = (namechainres)[10:11]
        resid= (namechainres)[15:]
        resid = int(resid)
        length = int(length)    
        os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/ids_xray")
        structure = PDBParser().get_structure(name,'pdb'+name+".ent")
        curves = []
        for residue_middle in range(resid+2, resid+length-2):
            c_vectors = []
            for residue in range(residue_middle-2,residue_middle+3,2):
                try:
                        c_alpha = structure[0][chain][residue]['CA']
                        c_alpha_angle = c_alpha.get_vector()
                        print(c_alpha.get_parent().get_resname(),residue,c_alpha_angle)
                        c_vectors.append(c_alpha_angle)

                except KeyError:
                    return np.NaN  
            angle = calc_angle(c_vectors[0],c_vectors[1],c_vectors[2])
            print(degrees(angle))
            curves.append(degrees(angle))
        return np.mean(curves)


def bigcurve():
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    phipsidf = pd.read_csv("phi_xray_chains_sorted_rev4.csv",usecols = ['line_number','name_chain_res','length'])
    #phipsidf = phipsidf.sort_values(by ='length', ascending = False)
    phipsidf['curvature'] = phipsidf.apply(lambda row: calccurve(row['name_chain_res'],row['length']) if row['length'] > 4 else np.NaN, axis = 1)
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    phipsidf.to_csv('phi_xray_chains_sorted_curvature_rev4.csv')
    return phipsidf

#print(calccurve_general("4bz4:ChainA:TYR208",11))
print(bigcurve())
