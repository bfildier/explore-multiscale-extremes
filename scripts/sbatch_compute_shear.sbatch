#!/bin/bash

# Compute distribution of rainfall at each time index. Use SLURM_ARRAY_TASK_ID as time index, and loop over indices in relation table DYAMOND-SAM/TOOCAN segmentation mask

# ONLY WORKS IF all.nc not computed

#SBATCH --partition=zen16
#SBATCH -o logs/run_%x_%j_A%A_%a.o
#SBATCH -e logs/run_%x_%j_A%A_%a.e
#SBATCH --nodes=1
#SBATCH --time=06:00:00
#SBATCH --job-name=lowshear
#SBATCH --mem=10g
#SBATCH --ntasks-per-node=1

# job information 
cat << EOF
------------------------------------------------------------------
Job submit on $SLURM_SUBMIT_HOST by $SLURM_JOB_USER
JobID=$SLURM_JOBID Running_Node=$SLURM_NODELIST 
Node=$SLURM_JOB_NUM_NODES Task=$SLURM_NTASKS
------------------------------------------------------------------
EOF
# Begin of section with executable commands

echo 'source own .bashrc file'
source /home/bfildier/.bashrc

echo 'load conda environment'
conda activate pyLMD

SCRIPT_DIR=/home/bfildier/analyses/explore-multiscale-extremes/scripts/

# compute distribution at each index
srun -n 1 -c 1 --exclusive ${SCRIPT_DIR}/computeShear.sh

exit 0
