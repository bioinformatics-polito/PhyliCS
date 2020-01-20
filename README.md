# PhyliCS (TODO)
A pipeline for multi-sample copy-number variation (CNV) analysis on single-cell DNA sequencing data and intra-tumor heterogeneity quantification. 

## Installation and setup instructions

Download the .zip archive of the project and uncompress it or clone it via git, typing:

```git clone https://github.com/bioinformatics-polito/PhyliCS.git```

You need to create and activate a conda environment which will contain all the packages required to run PhyliCS. Then you can build the application by executing `build.sh`:
```
cd PhyliCS
ENV=<env_name> #substitute <env_name> with your preferred environment name
conda create --name $ENV --file phylics_env.txt
conda activate $ENV

./build.sh
```
Once you have built the application, add the directory containing the binaries to your `$PATH`. To do so, open `~/.bashrc` and add the following line:

```export PATH=<path/to/phylics>/PhyliCS/bin:$PATH```
