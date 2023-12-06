#!/bin/bash -e
#SBATCH --job-name=bh_polaron_var # job name (shows up in the queue)
#SBATCH --time=3-0:00:00           # Walltime (HH:MM:SS)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=0                    # memory/cpu (in MB)
#SBATCH --partition=short          # specify a partition
#SBATCH --hint=nomultithread       # don't use hyperthreading
#SBATCH --error=%A_%a.err
#SBATCH --output=%A_%a.out
#SBATCH --mail-user=matijacufar@gmail.com
#SBATCH --mail-type=END,FAIL
#SBATCH --array=1-50
t=0.$(printf %02d ${SLURM_ARRAY_TASK_ID})

srun julia -tauto variational.jl -t${t} $*
