#!/bin/bash

# Compute distribution of rainfall at each time index. Use SLURM_ARRAY_TASK_ID as time index, and loop over indices in relation table DYAMOND-SAM/TOOCAN segmentation mask

#SBATCH --partition=zen16
#SBATCH -o logs/run_%x_%j_A%A_%a.o
#SBATCH -e logs/run_%x_%j_A%A_%a.e
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --job-name=MCSage
#SBATCH --mem=10g
#SBATCH --ntasks-per-node=10
#SBATCH --array=832-1917:100%10
##SBATCH --ntasks-per-node=1

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

SCRIPT_DIR=/home/bfildier/analyses/explore-multiscale-extremes/scripts

# ARGUMENTS
# r : region (default = tropics)
region='tropics'
# m : mask (default = all)
mask_ref='all'

# mask_xy='allXY'
mask_xy='ocean'
mask_T='allT'
# mask_T='min5hr'

while getopts r:m flag
do
        case "${flag}" in
                r) region=${OPTARG}
                         ;;
		m) mask=${OPTARG}
			 ;;
                *) echo "Invalid option: -$flag" ;;
        esac
done

echo mask_ref = "${mask_ref}"
echo mask_xy = "${mask_xy}"
echo mask_T = "${mask_T}"

# compute distribution at each index
srun -n 1 -c 1 --exclusive python ${SCRIPT_DIR}/computeMcsAgeOnPrDistributionEachT.py -i ${SLURM_ARRAY_TASK_ID} -n 100 -r ${region} -m ${mask_ref} --mask_xy ${mask_xy} --mask_T ${mask_T}
#srun -n 1 -c 1 --exclusive python ${SCRIPT_DIR}/computeMcsAgeOnPrDistributionEachT.py -i 901 -n 1 -r ${region} -m ${mask_ref} --mask_xy ${mask_xy}

exit 0
