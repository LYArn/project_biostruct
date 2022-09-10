import numpy as np
import pandas as pd

DOPE = pd.read_csv('dope.par.txt', ' ', header=None)

#Create PDBParser object
with open("1gcn.pdb", 'r') as f:
    atlist = {}
    atnum = []
    atname = []
    atcoordx = []
    atcoordy = []
    atcoordz = []
    atcoord = [] 
    for line in f:
        if line.startswith('ATOM'):
            if line[13:17].strip() == 'CA':
                atnum.append(line[6:12].strip())
                atname.append(line[17:20].strip())
                atcoordx.append(float(line[30:38]))
                atcoordy.append(float(line[38:46]))
                atcoordz.append(float(line[46:54]))


for i in range(len(atnum)):
    atcoord.append(np.array([atcoordx[i], atcoordy[i], atcoordz[i]]))
    atlist[atnum[i] + atname[i]] = [atcoord[i]]


#Calculate distance between CA atoms i and all atoms, starting from i+3
ca_dist = []
ca_dist_asso = {}
for i in range(len(atlist)-2):
        for n in range(3, len(atlist)-i):
            diff_coord = atcoord[i] - atcoord[i+n]
            ca_dist.append(np.sqrt(np.sum(diff_coord * diff_coord)))


# #Associate pair of CA with their distance
# for i in range(len(ca_dist)):
#     key_name = str(i+1) + ca_name[i] + '_' + str(i+3) + ca_name[i + 2]
#     ca_dist_asso[key_name] = ca_dist[i]


# ca_dope_asso = {}
# value_ang = np.arange(0, 15.5, 0.5)
# value_dope = DOPE.iloc[:, 4:]
# value_asso = {}
# for i in range(len(DOPE)):
#     for n in range(30):
#         value_asso[value_ang[n]] = 


# for i in range(len(DOPE)):
#     if DOPE(i, 1) & DOPE(i, 3) == 'CA':
#         ca_dope_asso[DOPE(i, 0) + '_' + DOPE(i, 2)] = 
