import pandas as pd
import numpy as np
from Bio.PDB.PDBParser import PDBParser

file = "1gcn.pdb"
file_name = "Gluc"

#Create PDBParser object
parser = PDBParser() 
prot = parser.get_structure(file_name, file)

#Obtain coordinates of CA atoms
for model in prot:
   for chain in model:
        ca_coord = []
        for residue in chain:
                ca_coord.append(residue["CA"].get_coord())

#Calculate distance between CA atoms i and i+2
ca_dist = []
for i in range(len(ca_coord)):
    if i >= (len(ca_coord)-2):
        break
    else: 
        diff_coord = ca_coord[i] - ca_coord[i+2]
        ca_dist.append(np.sqrt(np.sum(diff_coord * diff_coord)))



