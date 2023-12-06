#!/bin/bash -e
#SBATCH --time=14-0:00:00           # Walltime (HH:MM:SS)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=0                    # memory/cpu (in MB)
#SBATCH --partition=long          # specify a partition
#SBATCH --hint=nomultithread       # don't use hyperthreading
#SBATCH --error=%A_%a
#SBATCH --output=%A_%a
#SBATCH --mail-user=matijacufar@gmail.com
#SBATCH --mail-type=END,FAIL
#SBATCH --array=1-20
t=0.$(printf %02d $SLURM_ARRAY_TASK_ID)

srun julia -tauto bosons.jl -otscan_long --warmup=20000 --measure=500000 -N1e7 --chunk_size=1000 -t${t} $*
