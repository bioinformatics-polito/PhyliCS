#!/bin/bash

#check if conda env exists
#   if not, create it

ROOT="$(pwd)"
#CONDA_ROOT="$(conda info --base)"

#ENV="phylics"

# enable the 'conda' command from bash
#source $CONDA_ROOT/etc/profile.d/conda.sh

# Check if the environment exists
#conda activate $ENV
#if [ $? != 0 ]; then
#    conda create --name $ENV --file $ROOT/phylics_env.txt
#    conda activate $ENV
#fi

mkdir -p bin
mkdir -p src

cp phylics/local/src/*.py src
cd bin

for p in $ROOT/src/*.py; do
    NOPATH=$(basename -- "$p")
    NOEXT="${NOPATH%.*}"
    ln -s $p $NOEXT
done 

#compile ginkgo
cd $ROOT/phylics/local/src/ginkgo 
make 2> /dev/null

cd $ROOT/phylics/local/src/ginkgo/genomes/scripts
make 2> /dev/null


#scripts cli genomes
cp $ROOT/phylics/local/src/ginkgo/cli/ginkgo.sh $ROOT/bin/
cp -r $ROOT/phylics/local/src/ginkgo/scripts $ROOT/
cp -r $ROOT/phylics/local/src/ginkgo/cli $ROOT/
cp -r $ROOT/phylics/local/src/ginkgo/genomes $ROOT/
