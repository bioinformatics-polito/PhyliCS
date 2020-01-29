# PhyliCS
A pipeline for multi-sample copy-number variation (CNV) analysis on single-cell DNA sequencing data and intra-tumor heterogeneity quantification. 

## Requirements

- Platform: 
    - linux-64
- Software:
    - python >= 3.6
    - conda >= 4.3

## Installation and setup instructions

### From source code
Download the archive containing the release version of the project and uncompress it:

```
wget -O PhyliCS-1.0.0.tar.gz https://github.com/bioinformatics-polito/PhyliCS/archive/v1.0.1.tar.gz
tar xvzf PhyliCS-1.0.0.tar.gz
```

You need to create and activate a conda environment which will contain all the packages required to run PhyliCS. Then you can build the application by executing `build.sh`:
```
cd PhyliCS
export PHYLICS_ENV=<env_name> #substitute <env_name> with your preferred environment name
conda create --name $PHYLICS_ENV --file phylics_env.txt
conda activate $ENV

./build.sh
```
Once you have built the application, add the directory containing the binaries to your `$PATH`. To do so, open `~/.bashrc` and add the following line:

```export PATH=<path/to/phylics>/PhyliCS/bin:$PATH```

### From Bioconda

Add Bioconda channel:

```conda config --add channels bioconda```

Create an environment containing PhyliCS and all of its dependencies:

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

    cd data/navin_primary 
    wget https://github.com/bioinformatics-polito/PhyliCS/blob/master/data/navin_primary/SegCopy
    wget https://github.com/bioinformatics-polito/PhyliCS/blob/master/data/navin_primary/results.txt

    cd ../navin_metastasis
    wget https://github.com/bioinformatics-polito/PhyliCS/blob/master/data/navin_metastasis/SegCopy
    wget https://github.com/bioinformatics-polito/PhyliCS/blob/master/data/navin_metastasis/results.txt

    cd ../..

    wget https://raw.githubusercontent.com/bioinformatics-polito/PhyliCS/master/phylics.ipynb

    #be sure to have activated phylics conda environment
    jupyter notebook
```
