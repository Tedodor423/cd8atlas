#!/bin/bash

# Resources:
#SBATCH --time=0-02:00:00  # DAYS-HOURS:MINUTES:SECONDS
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --partition=short

# Environment:
#SBATCH --export=NONE

# Email when job completes or fails
#SBATCH --mail-type=END,FAIL

# What to run:
module load python-cbrg

# custom pip packages
export PYTHONUSERBASE=/ceph/project/sharmalab/dnimrich/my-python
export PYTHONPATH=/ceph/project/sharmalab/dnimrich/my-python/lib/python3.11/site-packages${PYTHONPATH:+:$PYTHONPATH}


python /ceph/project/sharmalab/dnimrich/cd8atlas/code/test_pipeline.py
