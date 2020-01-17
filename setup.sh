#!/bin/bash

#check if conda env exists
#   if not, create it

#if ! which command > /dev/null; then
#    echo -e "Command not found! Install? (y/n) \c"
#    read
#    if "$REPLY" = "y"; then
#        sudo apt-get install command
#    fi
#fi
ROOT="$(pwd)"
#SCTOOLS=$ROOT/SCTools
#PHILICS=$ROOT/phylics


ENV="phylics"

# Check if the environment exists
conda activate $ENV
if [ $? -neq 0 ]; then
    # Create the environment and activate
    echo "Conda env '$ENV' doesn't exist."
    conda create --name $ENV --file $ROOT/phylics_env.txt
    conda activate $ENV
fi
   

#if [ ! -d $SCTOOLS ] ; then
#    #download sctools and build it
#    git clone https://github.com/bioinformatics-polito/SCTools.git
#    cd SCTools
#    mkdir build
#    cd build
#    cmake .. -DCMAKE_BUILD_TYPE=Release
#    make

#    cd $ROOT

#fi

#export SCTOOLS_DEMUX=$SCTOOLS/build/apps/sctools_demultiplex

#compile ginkgo
cd $ROOT/phylics/local/src/ginkgo
make 
cd $ROOT
