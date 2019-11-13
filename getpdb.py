from Bio.PDB import *


def getpdbid(sourcetxt):
    pdb1 = PDBList()
    pdbidlist = [line.rstrip() for line in open(sourcetxt,"r")]
    totallen = pdbidlist.__len__()
    for i,pdbid in enumerate(pdbidlist):
        cutpdbid = pdbid[0:4]
        print(cutpdbid)
        print("File %s out of %s" % (i+1,totallen))
        pdb1.retrieve_pdb_file(cutpdbid, pdir = "ids_xray",file_format = "pdb", overwrite = False)	

sourcedir = "C:/Users/Shaheer Rizwan/Documents/ramachandran/ids_xray.txt"
getpdbid(sourcedir)

