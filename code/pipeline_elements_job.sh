#!/bin/bash

# Resources:
#SBATCH --time=0-02:00:00  # DAYS-HOURS:MINUTES:SECONDS
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --partition=short

# Environment:
#SBATCH --export=NONE

# Email when job completes or fails
#SBATCH --mail-type=END,FAIL

# What to run:
module load python-cbrg
module load R-cbrg

# custom pip packages
export PYTHONUSERBASE=/ceph/project/sharmalab/dnimrich/my-python
export PYTHONPATH=/ceph/project/sharmalab/dnimrich/my-python/lib/python3.11/site-packages${PYTHONPATH:+:$PYTHONPATH}


python /ceph/project/sharmalab/dnimrich/cd8atlas/code/pipeline_elements.py
