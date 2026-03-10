#!/bin/bash

# Resources:
#SBATCH --time=0-00:10:00  # DAYS-HOURS:MINUTES:SECONDS
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --partition=test

# Environment:
#SBATCH --export=NONE

# Email when job completes or fails
#SBATCH --mail-type=END,FAIL

# What to run:
module load python-cbrg

# custom pip packages
export PYTHONUSERBASE=/ceph/project/sharmalab/dnimrich/my-python


python /ceph/project/sharmalab/dnimrich/cd8atlas/test_pipeline.py