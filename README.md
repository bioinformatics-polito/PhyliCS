# PhyliCS
A Python library to explore scCNA data and quantify spatial tumor heterogeneity

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
Download the built distribution file. 
You can do it manually by clicking [here](https://github.com/bioinformatics-polito/PhyliCS/raw/master/dist/phylics-1.0.0-py3-none-any.whl).

Install PhyliCS by typing:
```
pip install /path/to/phylics-1.0.0-py3-none-any.whl
```

### From PyPI

Installation trough PyPI is not available yet. We will upload PhyliCS to PyPI only after paper acceptance.

## Usage


## Case studies
Results discussed in the paper are all stored in a dedicated [repository](https://github.com/bioinformatics-polito/PhyliCS_usage) and summarized by means of a pair of  jupyter notebooks accessible trhough:
- [breast tumor](https://github.com/bioinformatics-polito/PhyliCS_usage/blob/main/breast_tumor.ipynb) 
- [lung tumor](https://github.com/bioinformatics-polito/PhyliCS_usage/blob/main/lung_tumor.ipynb)