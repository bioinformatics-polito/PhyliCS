# PhyliCS
a Python library to explore scCNA data and quantify spatial tumor heterogeneity

## Requirements

- python >=3.8, <3.9
- numpy>=1.19.5
- scipy>=1.6.0
- pandas>=1.1.3 
- matplotlib>=3.3.1
- seaborn>=0.11.1 
- scikit-learn>=0.24
- statsmodels>=0.12.0
- anndata>=0.7.5 
- typing
- umap-learn>=0.4.6
- IPython>=7.19.0
- hdbscan>=0.8.26
- joblib>=1.0.0

## Installation and setup instructions

If you do not have a working installation of Python 3.8.X (not 3.9), consider installing Miniconda (see [Miniconda](https://docs.conda.io/en/latest/miniconda.html) website). Then run:

```
conda create --name py38 python=3.8
conda activate py38
```

### From dist file
Download the built distribution file. You can do it manually by navigating to PhyliCS/dist/phylics-1.0.0-py3-none-any.whl and clicking on 'Download'.

Alternatively, you can do it from the command line by typing:

```
wget https://github.com/bioinformatics-polito/PhyliCS/raw/master/dist/phylics-1.0.0-py3-none-any.whl
```
Install PhyliCS by typing:
```
pip install /path/to/phylics-1.0.0-py3-none-any.whl
```

### From PyPI

Installation trough PyPI is not available yet. We will upload PhyliCS to PyPI only after paper acceptance.

## Test
Here we present the code to reproduce the results we have presented in our paper.

### Download the dataset

First you need to download the input data. We have provided the CNV calls (SegCopy files) and a statistics file (results.txt) for each dataset.

Go to the location where you want to store the data and type:

```
mkdir -p data/breast/breastA data/breast/breastB data/breast/breastC data/breast/breastD data/breast/breastE 
mkdir -p data/lung/primary data/lung/metastasis
    
cd data/breast/breastA
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/data/breast/breastA/SegCopy
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/data/breast/breastA/results.txt
    
cd ../breastB
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/data/breast/breastB/SegCopy
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/data/breast/breastB/results.txt
    
cd ../breastC
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/data/breast/breastC/SegCopy
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/data/breast/breastC/results.txt

cd ../breastD
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/data/breast/breastD/SegCopy
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/data/breast/breastD/results.txt

cd ../breastE
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/data/breast/breastE/SegCopy
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/data/breast/breastE/results.txt

cd ../../..

cd data/lung/primary 
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/data/lung/primary/SegCopy
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/data/lung/primary/results.txt

cd ../metastasis
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/data/lung/metastasis/SegCopy
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/data/lung/metastasis/results.txt

cd ../../..
    
```
### Run PhyliCS
