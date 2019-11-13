import math
from Bio.PDB import *
import os
import pandas as pd
import numpy as np


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

def ramachandran_type(residue, next_residue) :
    """Expects Bio.PDB residues, returns ramachandran 'type'

    If this is the last residue in a polypeptide, use None
    for next_residue.

    Return value is a string: "General", "Glycine", "Proline"
    or "Pre-Pro".
    """
    if residue.resname.upper()=="GLY" :
        return "Glycine"
    elif residue.resname.upper()=="PRO" :
        return "Proline"
    elif next_residue is not None \
    and next_residue.resname.upper()=="PRO" :
        #exlcudes those that are Pro or Gly
        return "Pre-Pro"
    else :
        return "General"


def get_chains():
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    df = pd.read_csv("ids_xray.txt",header = None, names=['chain'])
    df.set_index(df.chain.str.lower().str[0:4],inplace = True)
   
    df['chain'] = df.chain.str[4:]
    df.index.names = ['ids']
    #df2 = df.groupby(level = 'ids').agg({'chain':','.join}) - gives comma separated chain letters in each column
    df2 = df.groupby(level = 'ids')['chain'].agg({'chain': lambda x: list(x)})
    df.index.name = None
    return df2
   
  
    



print("loading Bio.PDB and the PDB file...")
os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/ids_xray")
df_chain = get_chains()

def unused():
    #root = "C:/Users/Shaheer Rizwan/Documents/ramachandran/ids_xray"
    #for path, subdirs, files in os.walk(root):
    #    for name in files:
    #        os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/ids_xray")
    #        structure = PDBParser().get_structure(name[:-4],name)
    #        print("About to save angles to file...")       
    #        for chains in df_chain.loc[[name[3:7],'chain']]:
    #            os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/phipsi_xray_chains")
    #            output_file = open("%s_chain%s_biopython.tsv" % (name[:-4], chains),"w")
    #            chain = structure[0][chains]
    #            print ("Chain %s" % str(chain.id))
    #            polypeptides = CaPPBuilder().build_peptides(chain)
    #            for poly_index, poly in enumerate(polypeptides) :
    #                phi_psi = poly.get_phi_psi_list()
    #                for res_index, residue in enumerate(poly) :
    #                    phi, psi = phi_psi[res_index]
    #                    if phi and psi :
    #						#Don't write output when missing an angle
    #                        output_file.write("%s
    #                        :Chain%s:%s%i\t%f\t%f\t%s\n" \
    #                        % (name[:-4], str(chain.id), residue.resname,
    #                        residue.id[1], degrees(phi), degrees(psi),
    #                        ramachandran_type(residue, poly[res_index+1])))
    #            output_file.close()
    #            print("Done")
    #print("DIHEDRALS FOR ALL PDB FILES HAVE BEEN CALCULATED yay")
    pass


def get_dihedrals(id_list):    
    for name,chain_list in id_list.iterrows():
        try:
            os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/ids_xray")
            pdbid = 'pdb'+name+'.ent'
            #pdbid = 'pdb4ybb.ent'
            structure = PDBParser().get_structure(name,pdbid)
            print("About to save angles to file...")
            for chains in chain_list:
                for chainid in chains:
                    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/phipsi_xray_chains")
                    output_file = open("pdb%s_chain%s_biopython.tsv" % (name, chainid),"w")
                    chain = structure[0][chainid]
                    print ("Chain %s" % str(chain.id))
                    polypeptides = CaPPBuilder().build_peptides(chain)
                    for poly_index, poly in enumerate(polypeptides) :
                        phi_psi = poly.get_phi_psi_list()
                        for res_index, residue in enumerate(poly) :
                            phi, psi = phi_psi[res_index]
                            if phi and psi :
						        #Don't write output when missing an angle
                                output_file.write("%s:Chain%s:%s%i\t%f\t%f\t%s\n" \
                                % (name, str(chain.id), residue.resname,
                                residue.id[1], degrees(phi), degrees(psi),
                                ramachandran_type(residue, poly[res_index+1])))
                    output_file.close()
                    print("Done")
        except FileNotFoundError:
            continue
        
get_dihedrals(df_chain)
