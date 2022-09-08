import numpy as np
import pandas as pd
from Bio.PDB.PDBParser import PDBParser

file = "1gcn.pdb"
file_name = "Gluc"
DOPE = pd.read_csv('dope.par.txt', ' ', header=None)

#Create PDBParser object
parser = PDBParser() 
prot = parser.get_structure(file_name, file)

#Obtain residues' name and CA atoms' coordinates
for model in prot:
   for chain in model:
        ca_coord = []
        ca_name = []
        for residue in chain:
            if residue.has_id('CA'):
                ca_coord.append(residue['CA'].get_coord())
                ca_name.append(residue.get_resname())

#Calculate distance between CA atoms i and i+2
ca_dist = []
ca_dist_asso = {}
for i in range(len(ca_coord)):
    if i >= (len(ca_coord)-2):
        break
    else: 
        diff_coord = ca_coord[i] - ca_coord[i+2]
        ca_dist.append(np.sqrt(np.sum(diff_coord * diff_coord)))

#Associate pair of CA with their distance
for i in range(len(ca_dist)):
    key_name = str(i+1) + ca_name[i] + '_' + str(i+3) + ca_name[i + 2]
    ca_dist_asso[key_name] = ca_dist[i]


ca_dope_asso = {}
value_ang = np.arange(0, 15.5, 0.5)
value_dope = DOPE.iloc[:, 4:]
value_asso = {}
for i in range(len(DOPE)):
    for n in range(30):
        value_asso[value_ang[n]] = 


# for i in range(len(DOPE)):
#     if DOPE(i, 1) & DOPE(i, 3) == 'CA':
#         ca_dope_asso[DOPE(i, 0) + '_' + DOPE(i, 2)] = 
