import numpy as np
import pandas as pd
import random
import copy
from math import sqrt


class Structure:
    def __init__(self):
        self.atname = []
        self.atcoord = [] 

    def get_coord(self, file):
        """ Parsing of .pdb file
            Return residues CA's coordinates and names"""
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
    value_ang = np.arange(0.5, 15.5, 0.5)
    value_dope = {}

    def __init__(self):
        Structure.__init__(self)
        self.ca_dist = []
        self.ca_name = []
        self.result = {}

    def init_dope(self, DOPE):
        #Create dataframe of DOPE value with Angstrom distance as column's name
        col = ['res1', 'res2']
        for i in range(len(self.value_ang)):
            col.append(self.value_ang[i])
        self.value_dope = DOPE[(DOPE[1]=='CA') & (DOPE[3]=='CA')]
        self.value_dope = self.value_dope.drop([1, 3], axis = 1)
        self.value_dope = self.value_dope.reset_index(drop=True)
        self.value_dope = self.value_dope.set_axis(col, axis = 1)
    
    def distance(self):
        """ Calculate distance between CA residue i and other residues i, with i > 3
            Return list of distances between CA and list of CA pairs (names with residues' number)"""
        for i in range(len(self.atname)-2):
            for n in range(3, len(self.atname)-i):
                diff_coord = self.atcoord[i] - self.atcoord[i+n]
                self.ca_dist.append(np.sqrt(np.sum(diff_coord * diff_coord)))
                self.ca_name.append(self.atname[i] + '_' + self.atname[i+n])

    def potential_score(self):
        self.asso_name_index = {}
        self.asso_name_dist = {}
        self.asso_name_ang = {}

        self.init_dope(DOPE)
        #Select row of interest in value_dope for each experimental distance
        for i in range(len(self.ca_name)):
            for index, row in self.value_dope.iterrows():
                if self.ca_name[i] == f"{row['res1']}_{row['res2']}":
                    self.asso_name_index[f"{i+1}{self.ca_name[i]}"] = index
        
        #Distance associated with the pair names
        for i in range(len(self.ca_dist)):
            self.asso_name_dist[f'{i+1}{self.ca_name[i]}'] = self.ca_dist[i]
        
        #Select correct distance in Angstrom for each experimental distance 
        for n in self.asso_name_index.keys():
            for i in range(len(self.value_ang)):
                if self.asso_name_dist[n] > self.value_ang[-1]:
                    self.asso_name_ang[n] = self.value_ang[-1]
                    break
                elif self.asso_name_dist[n] > self.value_ang[i]:
                    pass
                elif self.asso_name_dist[n] < self.value_ang[i]:
                    self.asso_name_ang[n] = self.value_ang[i]
                    break
                
        #Values of index and col being determined, associating pdf with CA pairs
        for i in self.asso_name_index.keys():
            self.result[i] = self.value_dope.loc[self.value_dope.index[self.asso_name_index[i]], self.asso_name_ang[i]]

def randomize_value(rand, rand_total):
    #Randomize list of CA to obtain a ref value to compare
    temp = list(zip(rand.atname, rand.atcoord))
    random.shuffle(temp)
    rand.atname, rand.atcoord = zip(*temp)
    rand.atname, rand.atcoord = list(rand.atname), list(rand.atcoord)
    rand.distance()
    rand.potential_score()
    rand_total.append(sum(rand.result.values()))



if __name__ == "__main__":
    DOPE = pd.read_csv('data/dope.par.txt', ' ', header=None)
    print("Input .pdb file's name :")

    with open(f'data/1gcn.pdb', 'r') as file:
        prot = Residue() #Create instance of class Residues
        prot.get_coord(file) #Parse glucagon .pdb file
    
        prot.distance() #Calculate distance between each CA pair
        prot.potential_score() #Associate pdf value with CA pair distances
        prot_total = sum(prot.result.values())
        
        rand = copy.deepcopy(prot)
        rand_total = []
        times = 10
        for i in range(times):
            randomize_value(rand, rand_total)

    #Z-score
    mean = sum(rand_total) / len(rand_total)
    print(rand_total)
    diff = [(value - mean)**2 for value in rand_total]
    sum_of_diff = sum(diff)
    standard_deviation = (sum_of_diff / (len(rand_total) - 1 )) ** 0.5
    print(mean, standard_deviation)

    #Statistical comparison
    score = (mean - prot_total)/standard_deviation
    print(score)

    # #Save results in a txt file
    with open('results/compil_results.txt', 'w') as result_file:
        result_file.write(f"File name : {f}\n")
        result_file.write(f"Sum of DOPE values for original sequence = {prot_total}\n")
        result_file.write(f"Mean and standard deviation calculated from randomized sequence ({times}x) = {mean}, {standard_deviation}\n")
        result_file.write(f"The score of the original sequence compared to randomized sequences : {score}\n")
        print('compil_results updated !')
