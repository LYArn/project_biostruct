import numpy as np
import pandas as pd

class Structure:
    def __init__(self):
        self.atname = []
        self.atcoord = [] 

    def pars(self, file):
        """ Parsing of .pdb file
            Return CA residues' coordinates and names"""
        x = []
        y = []
        z = []

        for line in file:
            if line.startswith('ATOM') and line[13:17].strip() == 'CA':
                self.atname.append(line[17:20].strip())
                x.append(float(line[30:38]))
                y.append(float(line[38:46]))
                z.append(float(line[46:54]))

        for i in range(len(self.atname)):
            self.atcoord.append(np.array([x[i], y[i], z[i]]))

class Residue(Structure):
    value_dope = {}

    def __init__(self):
        Structure.__init__(self)
        self.ca_dist = []
        self.ca_name = []
        self.result = {}


    def distance(self):
        """ Calculate distance between CA residue i and other residues i, with i > 3
            Return list of distances between CA and list of CA pairs (names with residues' number)"""
        for i in range(len(self.atname)-2):
            for n in range(3, len(self.atname)-i):
                diff_coord = self.atcoord[i] - self.atcoord[i+n]
                self.ca_dist.append(np.sqrt(np.sum(diff_coord * diff_coord)))
                self.ca_name.append(self.atname[i] + '_' + self.atname[i+n])

    def cal_pdf(self,DOPE):
        self.distance()
        value_ang = np.arange(0.5, 15.5, 0.5)
        self.asso_name_index = {}
        self.asso_name_dist = {}
        self.asso_name_ang = {}

        #Create dataframe of DOPE value with Angstrom distance as column's name
        col = ['res1', 'res2']
        for i in range(len(value_ang)):
            col.append(value_ang[i])
        value_dope = DOPE[(DOPE[1]=='CA') & (DOPE[3]=='CA')]
        value_dope = value_dope.drop([1, 3], axis = 1)
        value_dope = value_dope.reset_index(drop=True)
        value_dope = value_dope.set_axis(col, axis = 1)

        #Select row of interest in value_dope for each experimental distance
        for i in range(len(self.ca_name)):
            for index, row in value_dope.iterrows():
                if self.ca_name[i] == f"{row['res1']}_{row['res2']}":
                    self.asso_name_index[f"{i+1}{self.ca_name[i]}"] = index
        
        #Distance associated with the pair names
        for i in range(len(self.ca_dist)):
            self.asso_name_dist[f'{i+1}{self.ca_name[i]}'] = self.ca_dist[i]
        
        #Select correct distance in Angstrom for each experimental distance 
        for n in self.asso_name_index.keys():
            for i in range(len(value_ang)):
                if self.asso_name_dist[n] > value_ang[-1]:
                    self.asso_name_ang[n] = value_ang[-1]
                    break
                elif self.asso_name_dist[n] > value_ang[i]:
                    pass
                elif self.asso_name_dist[n] < value_ang[i]:
                    self.asso_name_ang[n] = value_ang[i]
                    break
                
        #Values of index and col being determined, associating pdf with CA pairs
        for i in self.asso_name_index.keys():
            self.result[i] = value_dope.loc[value_dope.index[self.asso_name_index[i]], self.asso_name_ang[i]]

if __name__ == "__main__":
    DOPE = pd.read_csv('dope.par.txt', ' ', header=None)
    with open("1gcn.pdb", 'r') as file:
        gluc = Residue() #Create instance of class Residues
        gluc.pars(file) #Parse glucagon .pdb file
        gluc.cal_pdf(DOPE) #Associate pdf value with CA pair distances
        gluc_total = sum(gluc.result.values())
        print(f'Total energy score is {gluc_total}')
