#!/bin/bash



## Basic Run info 

#SBATCH --job-name=${JOB_NAME}

#SBATCH --partition=sandbox

#SBATCH --time=0-04:00:00 ## time format is DD-HH:MM:SS

#SBATCH --cpus-per-task=1

#SBATCH --mem=6400 ## max amount of memory per node you require



## Path to output and error file 

#SBATCH --error=${ERROR_FILE} ## %A - filled with jobid

#SBATCH --output=${OUTPUT_FILE} ## %A - filled with jobid



## Useful for remote notification

## make sure you update the bottom bits with your email before uncommenting 

#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE

#SBATCH --mail-user=kp779@hawaii.edu



## initializes the fenicsproject env

source ~/.zshrc

conda activate optoproject



## Basic run with defaults:

## Can run with sbatch level_set_op.slurm

##python3 mechanism_kp.py inverter

cd /home/kp779/BRIDGE_Connection/src/BRIDGE_Connection 

python3 compliance.py cantilever

##Run FEniCS simulations with custom parameters.


