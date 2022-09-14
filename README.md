# Assessment of the quality of protein's 3D model
## General information

The script of this project use a PDB extension file and a DataFrame containing Discrete Optimized Protein Energy (DOPE), statistical potential values. These values are precalculated and are downloaded from this [link](www.dsimb.inserm.fr/~gelly/data/dope.par).

By inputting a PDB file, we calculate the sum of DOPE values related to the interaction of each amino acid between each other. We assume these interactions are only between their alpha carbon.


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
conda activate dope_project-env
```

## Using script
