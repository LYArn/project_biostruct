import numpy as np
import pandas as pd
import random
import copy
import argparse
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
        """Create dataframe of DOPE value with Angstrom distance as column's name"""
        col = ['res1', 'res2']
        for i in range(len(self.value_ang)):
            col.append(self.value_ang[i])
        self.value_dope = DOPE[(DOPE[1]=='CA') & (DOPE[3]=='CA')]
        self.value_dope = self.value_dope.drop([1, 3], axis = 1)
        self.value_dope = self.value_dope.reset_index(drop=True)
        self.value_dope = self.value_dope.set_axis(col, axis = 1)
    
    def distance(self):
        """ Calculate distance between CA residue i and other residues i, with i > 3
            Return list of distances between CA and list of CA pairs 'names with residues' number
            Used in Method potential_score()"""
        for i in range(len(self.atname)-2):
            for n in range(3, len(self.atname)-i):
                diff_coord = self.atcoord[i] - self.atcoord[i+n]
                self.ca_dist.append(np.sqrt(np.sum(diff_coord * diff_coord)))
                self.ca_name.append(self.atname[i] + '_' + self.atname[i+n])

    def potential_score(self):
        """Using Dataframe of DOPE value from init_dope() and distances/names data from distance(), generate :
                A dictionary with CA pairs associated with row index (pairs' name) from the Dataframe
                A dictionary with CA pairs associated with columns values (Angstrom distances)
                Last dictionary using both dictionaries to associate CA pairs with the correct DOPE value"""
        self.asso_name_index = {}
        self.asso_name_dist = {}
        self.asso_name_ang = {}

        self.distance()
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
                
        #Values of index and col being determined, associating potential values with CA pairs
        for i in self.asso_name_index.keys():
            self.result[i] = self.value_dope.loc[self.value_dope.index[self.asso_name_index[i]], self.asso_name_ang[i]]

def randomize_value(rand, rand_total):
    """ Create values to compare DOPE value from the source sequence
        Using amino acids data from the source file, creation of a shuffled sequence.
        Output : list containing the sum of DOPE values from the randomized sequence"""
    rand.result.clear()
    rand.ca_dist.clear()
    rand.ca_name.clear()
    temp = list(zip(rand.atname, rand.atcoord))
    random.shuffle(temp)
    rand.atname, rand.atcoord = zip(*temp)
    rand.atname, rand.atcoord = list(rand.atname), list(rand.atcoord)
    rand.potential_score()
    rand_total.append(sum(rand.result.values()))

if __name__ == "__main__":
    DOPE = pd.read_csv('data/dope.par.txt', ' ', header=None)
    
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    args = parser.parse_args()

    with open(f'data/{args.filename}', 'r') as file:
        prot = Residue() #Create instance of class Residues
        prot.get_coord(file)
        prot.potential_score() #Associate DOPE value with CA pair distances
        prot_total = sum(prot.result.values())
        
        rand = copy.deepcopy(prot) #Copy data from the instance to rand without linking variables
        rand_total = []
        times = 10
        for i in range(times): #Looping 10 times to obtain 10 sum of DOPE values
            randomize_value(rand, rand_total)

    #Z-score calcul
    mean = sum(rand_total) / len(rand_total)
    diff = [(value - mean)**2 for value in rand_total]
    sum_of_diff = sum(diff)
    standard_deviation = (sum_of_diff / (len(rand_total) - 1 )) ** 0.5

    #Statistical comparison between experimental value (source sequence) and reference values (shuffled sequences)
    score = (mean - prot_total)/standard_deviation

    #Save results in a txt file
    with open('results/compil_results.txt', 'w') as result_file:
        result_file.write(f"File name : {args.filename}\n")
        result_file.write(f"Number of CA in this sequence : {len(prot.atname)}\n")
        result_file.write(f"Sum of DOPE values for original sequence = {round(prot_total, 2)}\n")
        result_file.write(f"List of ref values from shuffled sequences : {rand_total}\n")
        result_file.write(f"Mean and standard deviation calculated from randomized sequence ({times}x) = {round(mean, 2)}, {round(standard_deviation, 2)}\n")
        result_file.write(f"The score of the original sequence compared to randomized sequences : {round(score, 2)}\n")
        print('compil_results updated !')
