#!/bin/bash
# this is intended for running DPS jobs; the input directory is where a single file has been pulled because download=TRUE in the algorithm_config.yaml file

# This installs the python libs needed to run the script at the bottom

set -x

source activate python

unset PROJ_LIB

mkdir output

basedir=$( cd "$(dirname "$0")" ; pwd -P ) 

## Hard coded args for each run (if any; usually just output dir)

# Work dir is always from where your script is called
# Base dir is always the relative dir within the run*.sh script

# Absolute path here
# This PWD is wherever the job is run (where the .sh is called from) 
OUTPUTDIR="${PWD}/output"

python ${basedir}/../../lib/download_gee_tiles.py \
--tile ${1} \
--ys ${2} \
--ye ${3} \
--ds ${4} \
--de ${5} \
--output ${OUTPUTDIR}
