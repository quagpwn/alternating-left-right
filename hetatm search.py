import numpy as np
from Bio.PDB import *
import os
import pandas as pd
import ast
import time

print("loading Bio.PDB and the PDB file...")

def hetcalc(namechainres,length,type):
    name = (namechainres)[0:4]
    chain = (namechainres)[10:11]
    print(namechainres)
    resid= (namechainres)[15:]
    resid = int(resid)
    #print(resid)
    length = int(length)    
    hetlist =[]
    hetlist_simple = []
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/ids_xray")
    structure = PDBParser().get_structure(name,'pdb'+name+".ent")
    for residue in range(resid, resid+length):
        try:
            nitrogen = structure[0][chain][residue][type]
            #print(nitrogen.get_parent().get_resname(),residue)
            for hetres in structure[0][chain]:
                if hetres.id[0] != ' ':
                    for hetatm in hetres:
                        #print(hetatm.get_name())
                        dist = (nitrogen - hetatm)
                        if dist < 5.0 and hetatm.get_parent().get_resname() != "HOH" and hetatm.get_parent().get_resname() != "MSE" and hetatm.get_parent().get_resname() != "EDO" and hetatm.get_parent().get_resname() != "UNX":
                            hetlist.append([nitrogen.get_parent().get_resname()+str(residue),hetatm.get_name(),hetatm.get_parent().get_resname(), hetatm.get_parent().id[1], str(hetatm.get_serial_number()), str(dist)])
                            if not any(hetatm.get_parent().id[1] in sublist for sublist in hetlist_simple):
                                hetlist_simple.append([hetatm.get_parent().get_resname(),hetatm.get_parent().id[1]])
                            #print(dist,hetatm.get_name(),hetatm.get_parent().get_resname(), hetatm.get_serial_number(),"\n")
        except KeyError:
            continue
    return hetlist_simple,hetlist

#hetcalc_simple not to be used anymore, function has been integrated into hetcalc
def hetcalc_simple(namechainres,length,type): #gives just heteroatom molecules, no data about residues or atoms or serial number
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
            nitrogen = structure[0][chain][residue][type]
            #print(nitrogen.get_parent().get_resname(),residue)
            for hetres in structure[0][chain]:
                if hetres.id[0] != ' ':
                    for hetatm in hetres:
                        #print(hetatm.get_name())
                        dist = (nitrogen - hetatm)
                        if dist < 5.0 and hetatm.get_parent().get_resname() != "HOH" and hetatm.get_parent().get_resname() != "MSE" and hetatm.get_parent().get_resname() != "EDO":
                            if hetatm.get_parent().get_resname() not in hetlist:
                                hetlist.append([hetatm.get_parent().get_resname(),hetatm.get_parent().id[1]])
                            #print(dist,hetatm.get_name(),hetatm.get_parent().get_resname(), hetatm.get_serial_number(),"\n")
        except KeyError:
            continue
    return hetlist_simple,hetlist



def hetcalc_multi(namechainres,length, hetatm_simple,hetatm):#multiple hets? - compare N and O distances and see which set is legitimate, to be used with preexisting hetatm listings
    hetatm_simple = ast.literal_eval(hetatm_simple)
    hetatm = ast.literal_eval(hetatm)
    #print(hetatm_simple[0])
    name = (namechainres)[0:4]
    chain = (namechainres)[10:11]
    #print(namechainres)
    resid = (namechainres)[15:]
    resid = int(resid)
    #print(resid)
    length = int(length)
    if len(hetatm_simple) > 1:
        n_list = hetatm
        o_list = hetcalc(namechainres,length,'O')[1]
        print(len(o_list),'skrt')
        if len(o_list) == 0:
            hetlist = []
            hetlist_simple = []
            for ls in n_list:
                if float(ls[5]) <= 3: #add any hets that are actually close
                   if not any(int(ls[3]) in sublist for sublist in hetlist_simple):
                                hetlist_simple.append([ls[2],int(ls[3])])
                   hetlist.append(ls+['N'])

            hetdist = {} #also add the single closest het if not yet added
            for ls in n_list:
                if int(ls[3]) not in hetdist: #if key is new, create list containing first distance
                    hetdist[int(ls[3])] = [float(ls[5])]
                else: #if key is not new, add another distance 
                    hetdist[int(ls[3])].append(float(ls[5]))
                    #print(hetdist[ls[3]])            
            mean_hetdist = {} #take all the distances of each resn and get the mean in a new dict
            for key in hetdist:
                mean_hetdist[key] = np.mean(hetdist[key])
            print(hetdist)
            print(mean_hetdist)
            het_resn = min(mean_hetdist, key=mean_hetdist.get) #get the resn number of the one with the lowest mean distance
            print(het_resn,type(het_resn))
            for ls in n_list:
                if het_resn == int(ls[3]):
                    if not any(int(ls[3]) in  subl for subl in hetlist_simple):
                        hetlist_simple.append([ls[2],int(ls[3])])
                        print('yeet')
                    hetlist.append(ls+['N'])

            return hetlist_simple,hetlist

        if len(o_list) > 0:
            #print(o_list)
            hetlist = []
            hetlist_simple = []
            for ls in n_list:
                if float(ls[5]) <= 3: #add any hets that are actually close in n_list
                    if not any(int(ls[3]) in sublist for sublist in hetlist_simple):
                                hetlist_simple.append([ls[2],int(ls[3])])
                    hetlist.append(ls+['N'])
            for ls in o_list:
                if float(ls[5]) <= 3: #add any hets that are actually close in o_list
                    if not any(int(ls[3]) in sublist for sublist in hetlist_simple):
                                hetlist_simple.append([ls[2],int(ls[3])])
                    hetlist.append(ls+['O'])

            n_hetdist = {} #get the n het distances
            for ls in n_list:
                if int(ls[3]) not in n_hetdist: #if key is new, create list containing first distance
                    n_hetdist[int(ls[3])] = [float(ls[5])]
                else: #if key is not new, add another distance 
                    n_hetdist[int(ls[3])].append(float(ls[5]))
                    #print(hetdist[ls[3]])            
            n_mean_hetdist = {} #take all the distances of each resn and get the mean in a new dict 
            for key in n_hetdist:
                n_mean_hetdist[key] = np.mean(n_hetdist[key])

            o_hetdist = {}  #get the o het distances
            for ls in n_list:
                if int(ls[3]) not in o_hetdist: #if key is new, create list containing first distance
                    o_hetdist[int(ls[3])] = [float(ls[5])]
                else: #if key is not new, add another distance 
                    o_hetdist[int(ls[3])].append(float(ls[5]))
                    #print(hetdist[ls[3]])            
            o_mean_hetdist = {} #take all the distances of each resn and get the mean in a new dict
            for key in o_hetdist:
                o_mean_hetdist[key] = np.mean(o_hetdist[key])

            o_het_resn = min(o_mean_hetdist, key = o_mean_hetdist.get) #get single closest resn for o and n
            n_het_resn = min(n_mean_hetdist, key=n_mean_hetdist.get)
            if n_mean_hetdist[n_het_resn] < o_mean_hetdist[o_het_resn]: #add resn with the smallest distance to hetatm and hetatm_simple
                for ls in n_list:
                    if n_het_resn == int(ls[3]):
                        if not any(int(ls[3]) in  subl for subl in hetlist_simple):
                            hetlist_simple.append([ls[2],int(ls[3])])
                            print('yeet')
                        hetlist.append(ls+['N'])
            else:
                for ls in o_list:
                    if o_het_resn == int(ls[3]):
                        if not any(int(ls[3]) in  subl for subl in hetlist_simple):
                            hetlist_simple.append([ls[2],int(ls[3])])
                            print('yeeto')
                        hetlist.append(ls+['O'])
            return hetlist_simple,hetlist
    else:
        return hetatm_simple,hetatm


def gethet():
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    phipsidf = pd.read_csv("phi_xray_chains_sorted_curve_bfac_is_pos_cluster_conf_hetatmsimple_mse_edo_unx_25clust_TRIG_rev4.csv",usecols = ['line_number','name_chain_res','length','bfac_within_stddev','curvature','is_pos','cluster','conf_seq','hetatm','hetatm_simple'])
    #phipsidf = phipsidf.replace(r'\\n',' ', regex=True) 
    phipsidf = phipsidf.replace(r'\\n',' ', regex=True) #strip
    #phipsidf['hetatm'] = np.nan #fill with nans
    #phipsidf['hetatm_simple'] = np.nan #fill with nans
    phipsidf['name_chain_res'] = phipsidf['name_chain_res'].str.strip(r'\\n') #strip
    #phipsidf['hetatm'] = phipsidf.astype('object') #prep column by making into objects
    #phipsidf['hetatm_simple'] = phipsidf.astype('object') #prep column by making into objects
    phipsidf[['hetatm_simple','hetatm']] = phipsidf.apply(lambda row: hetcalc_multi(row['name_chain_res'], row['length'],row['hetatm_simple'],row['hetatm']),axis = 1,result_type = 'expand')
    #phipsidf[['hetatm_simple','hetatm']] = phipsidf.apply(lambda row: hetcalc(row['name_chain_res'], row['length'],'N'),axis = 1,result_type = 'expand')
    #print(phipsidf['hetatm'])
    os.chdir("C:/Users/Shaheer Rizwan/Documents/ramachandran/")
    phipsidf.to_csv('phi_xray_chains_sorted_curve_bfac_is_pos_cluster_conf_hetatmmulti_mse_edo_unx_25clust_TRIG_rev4.csv')
    #phipsidf.to_hdf('phipsi_xray_chains_hetatmlisted_no_water.h5',key='df',mode='w')
    #print(phipsidf)

start = time.time()
#print(hetcalc("4j5r:ChainA:ASP37",5,'N')[0])
#print("\n",hetcalc_multi('4lp8:ChainA:ALA95',7,"[['  K', 402], ['  K', 403], ['  K', 404], ['  K', 405]]","[['ILE97', 'K', '  K', 402, '2253', '4.842594'], ['ILE97', 'K', '  K', 403, '2254', '4.377421'], ['GLY98', 'K', '  K', 403, '2254', '4.5210195'], ['GLY98', 'K', '  K', 404, '2255', '4.436939'], ['TYR99', 'K', '  K', 404, '2255', '4.3555136'], ['TYR99', 'K', '  K', 405, '2256', '4.6970725'], ['GLY100', 'K', '  K', 405, '2256', '4.2039285']]"))
gethet()
end = time.time()
print ('Run time is:'+str(end - start))   
