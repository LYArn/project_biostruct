# Assessment of the quality of protein's 3D model using DOPE theory
## General information

The script of this project use a PDB extension file and a DataFrame containing Discrete Optimized Protein Energy (DOPE), statistical potential values. These values are precalculated and are downloaded from this [link](http://www.dsimb.inserm.fr/~gelly/data/dope.par). The theory used in this project is from [(Shen MY, 2006)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2242414/).

By inputting a PDB file, we calculate the sum of DOPE values related to the interaction of each amino acid between each other. We assume these interactions are only between their alpha carbon. 

As results, we obtain a score by comparing the sum of DOPE values with the z-score of randomized sequences (shuffled original sequence). After comparison, the higher (or lower, in negative value) the score is, the more likely this 3D conformation is to exist *in vivo*.


## Setting environment
### Download repository
```bash
git clone https://github.com/LYArn/project_biostruct.git
cd project_biostruct
```

### Install environment
Install conda.

Install mamba.

Create conda environment for the project :
```bash
mamba env create -f bin/environment.yml
```

Load conda environment : 
```bash
conda activate dope-env
```

## Using script

Call script with :
```bash
python scripts/main.py
```

Input the name of a PDB file from the folder `data` :

*(Run ~30 sec)*
```bash
1gcn.pdb
```
Or :

*(Run ~70 sec)*
```bash
1pdc.pdb
```

**Results of the run will be written in `compil_results.txt` in `results` folder**
