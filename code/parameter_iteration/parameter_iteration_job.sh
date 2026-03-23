#!/bin/bash

# Resources:
#SBATCH --time=0-20:00:00  # DAYS-HOURS:MINUTES:SECONDS
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --partition=short

# Environment:
#SBATCH --export=NONE

# Email when job completes or fails
#SBATCH --mail-type=END,FAIL

# What to run:
module load python-cbrg

# custom pip packages
export PYTHONUSERBASE=/ceph/project/sharmalab/dnimrich/my-python
export MPLCONFIGDIR=/tmp/cd8atlas_mpl
export NUMBA_CACHE_DIR=/tmp/cd8atlas_numba
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"

cd /ceph/project/sharmalab/dnimrich
python3 -m cd8atlas.code.parameter_iteration.parameter_iteration
