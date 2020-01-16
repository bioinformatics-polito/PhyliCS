#!/bin/bash

#===============================================================================
# Ginkgo extension to reclust multiple sample 
# data already analyzed through the pipeline
#===============================================================================

DIR_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd ../ && pwd )"
DIR_GENOMES=${DIR_ROOT}/genomes
DIR_SCRIPTS=${DIR_ROOT}/scripts

# ------------------------------------------------------------------------------
# Parse user input
# ------------------------------------------------------------------------------

# Sample usage
usage=" ---> Sample usage: /ginkgo/cli/reclust.sh --genome hg19 --input dir/to/result/files/ --results results.txt --segfixed SegFixed --segbreaks SegBreaks --binning variable_500000_101_bowtie [--clustlinkage ward] [--clustdist euclidean]  [--facs facs.txt] [--maskbadbins] [--masksexchrs]"


# Required parameters
unset DIR_INPUT
unset GENOME
unset BINNING
unset RESULTS
unset SEGFIXED
unset SEGBREAKS

# Optional parameters
CLUSTERING_DISTANCE="euclidean"
CLUSTERING_LINKAGE="ward"
FILE_FACS=""
MASK_SEXCHRS=0

# Parse user input
while [[ $# -gt 0 ]]; #1
do
    case "${1}" in
        # Required parameters
	--input         ) DIR_INPUT="$2"                                        ; shift; ;;
        --genome        ) GENOME="$2"                                           ; shift; ;;
        --binning       ) BINNING="$2"                                          ; shift; ;;
        --results       ) RESULTS="$2"                                          ; shift; ;;
	--segfixed      ) SEGFIXED="$2"                                         ; shift; ;;
	--segbreaks     ) SEGBREAKS="$2"                                        ; shift; ;;
	# Optional parameters
        --clustdist     ) CLUSTERING_DISTANCE="$2"                              ; shift; ;;
        --clustlinkage  ) CLUSTERING_LINKAGE="$2"                               ; shift; ;;
        --facs          ) FILE_FACS="$2"                                        ; shift; ;;
        --masksexchrs   ) MASK_SEXCHRS=1                                        ;        ;;
        *               ) echo "Warning: ignoring parameter ${1}"               ;        ;;
    esac
    shift
done

# Validate required parameters
DIR_INPUT=${DIR_INPUT?$usage}
GENOME=${GENOME?$usage}
BINNING=${BINNING?$usage}
RESULTS=${RESULTS?$usage}
SEGFIXED=${SEGFIXED?$usage}
SEGBREAKS=${SEGBREAKS?$usage}

# ------------------------------------------------------------------------------
# -- Setup variables
# ------------------------------------------------------------------------------

statFile=status.xml

# Make sure input folder is valid
[[ ! -d "${DIR_INPUT}" ]] && echo "Error: folder ${DIR_INPUT} doesn't exist" && exit

# Make sure input files exist
[[ ! -e "${RESULTS}" ]] && echo "Error: folder ${RESULTS} doesn't exist" && exit
[[ ! -e "${SEGFIXED}" ]] && echo "Error: folder ${SEGFIXED} doesn't exist" && exit
[[ ! -e "${SEGBREAKS}" ]] && echo "Error: folder ${SEGBREAKS} doesn't exist" && exit

 
# Genomes directory
DIR_GENOME=${DIR_GENOMES}/${GENOME}/"original"

# FACS file
FACS=$([[ -e "${FILE_FACS}" ]] && echo 1 || echo 0)

if [ "${FACS}" == 0 ];
then
    FILE_FACS=${DIR_INPUT}/"ploidyDummy.txt"
    touch ${FILE_FACS}
else
    # 
    uuid=$(uuidgen)
    # In case upload file with \r instead of \n (Mac, Windows)
    tr '\r' '\n' < ${FILE_FACS} > ${DIR_INPUT}/tmp-${uuid}
    mv ${DIR_INPUT}/tmp-${uuid} ${FILE_FACS}
    # 
    sed "s/.bed//g" ${FILE_FACS} | sort -k1,1 | awk '{print $1"\t"$2}' > ${DIR_INPUT}/tmp-${uuid} 
    mv ${DIR_INPUT}/tmp-${uuid} ${FILE_FACS}
fi

# ------------------------------------------------------------------------------
# -- Recreate Clusters
# ------------------------------------------------------------------------------
${DIR_SCRIPTS}/reclust.R ${DIR_SCRIPTS} ${DIR_INPUT} ${DIR_GENOME} ${RESULTS} ${SEGFIXED} ${SEGBREAKS} $statFile ${BINNING} ${CLUSTERING_LINKAGE} ${CLUSTERING_DISTANCE} ${FACS} ${FILE_FACS} ${MASK_SEXCHRS}
