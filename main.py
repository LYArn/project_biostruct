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
        ca_name = []
        for residue in chain:
            if residue.has_id('CA'):
                ca_coord.append(residue['CA'].get_coord())
                ca_name.append(residue.get_resname())

#Calculate distance between CA atoms i and i+2
ca_dist = []
ca_asso = {}
for i in range(len(ca_coord)):
    if i >= (len(ca_coord)-2):
        break
    else: 
        diff_coord = ca_coord[i] - ca_coord[i+2]
        ca_dist.append(np.sqrt(np.sum(diff_coord * diff_coord)))

#Associate pair of CA with their distance
for i in range(len(ca_dist)):
    ca_asso[ca_name[i] + '_' + ca_name[i + 2]] = ca_dist[i]
