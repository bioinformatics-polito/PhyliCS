# PhyliCS
A pipeline for multi-sample copy-number variation (CNV) analysis on single-cell DNA sequencing data and intra-tumor heterogeneity quantification. 

## Requirements

- Platform: 
    - linux-64
- Software:
    - python >= 3.6
    - conda >= 4.8.2

## Installation and setup instructions

### From source code
Clone the git archive:

```
git clone https://github.com/bioinformatics-polito/PhyliCS.git
```

You need to create and activate a conda environment which will contain all the packages required to run PhyliCS. Then you can build the application by executing `build.sh`:
```
cd PhyliCS
conda create --name phylics --file phylics_env.txt
conda activate phylics

./build.sh
```
Once you have built the application, add the directory containing the binaries to your `$PATH`. To do so, open `~/.bashrc` and add the following line:

```export PATH=<path/to/phylics>/PhyliCS/bin:$PATH```

### From Bioconda

Add Bioconda channel:

```conda config --add channels bioconda```

Create an environment containing PhyliCS and all its dependencies:

```
conda create --name phylics phylics
conda activate phylics
```

## Test
We have provided a Jupyter-notebook showing the code to test our application. It must be executed after activating the conda environment. 

If you have installed PhyliCS from the source code, then, just do:

```
cd PhyliCS
jupyter notebook
```

If you have installed it from conda, instead, then you need download the input data and the notebook from this repo and run the test. To do so, `cd` to the directory where you want to perform the test and type:

```
mkdir data
mkdir data/navin_primary data/navin_metastasis
mkdir data/10x_breastA data/10x_breastB data/10x_breastC

cd data/navin_primary 
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/data/navin_primary/SegCopy
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/data/navin_primary/results.txt

cd ../navin_metastasis
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/data/navin_metastasis/SegCopy
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/data/navin_metastasis/results.txt

cd ../..
    
cd data/10x_breastA
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/data/10x_breastA/SegCopy
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/data/10x_breastA/results.txt
    
cd ../10x_breastB
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/data/10x_breastB/SegCopy
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/data/10x_breastB/results.txt
    
cd ../10x_breastC
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/data/10x_breastC/SegCopy
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/data/10x_breastC/results.txt

cd ../..
    
wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/phylics.ipynb

#be sure you have activated phylics conda environment
jupyter notebook
```
