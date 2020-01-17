#!/bin/bash

#check if conda env exists
#   if not, create it

ROOT="$(pwd)"
CONDA_ROOT="$(conda info --base)"

ENV="phylics"

# enable the 'conda' command from bash
source $CONDA_ROOT/etc/profile.d/conda.sh

# Check if the environment exists
conda activate $ENV
if [ $? -neq 0 ]; then
    conda create --name $ENV --file $ROOT/phylics_env.txt
    conda activate $ENV
fi

mkdir -p bin
mkdir -p src

cp phylics/local/src/*.py src
for p in src; do
    NOPATH=$(basename -- "$p")
    NOEXT="${NOPATH%.*}"
    ln -s $p bin/$NOEXT;

#compile ginkgo
cd phylics/local/src/ginkgo & make
cd phylics/local/src/ginkgo/genomes/scripts & make

cp phylics/local/src/ginkgo/cli/ginkgo.sh $ROOT/bin/
cp -r phylics/local/src/ginkgo/* $ROOT/


