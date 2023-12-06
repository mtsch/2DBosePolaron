#!/bin/bash -e
#SBATCH --job-name=bh_polaron_nscan # job name (shows up in the queue)
#SBATCH --time=3-0:00:00           # Walltime (HH:MM:SS)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=0                    # memory/cpu (in MB)
#SBATCH --partition=short          # specify a partition
#SBATCH --hint=nomultithread       # don't use hyperthreading
#SBATCH --error=%A_%a
#SBATCH --output=%A_%a
#SBATCH --mail-user=matijacufar@gmail.com
#SBATCH --mail-type=END,FAIL
#SBATCH --array=3-9

srun julia -tauto bosons.jl -onscan_${SLURM_ARRAY_TASK_ID} --warmup=10000 --chunk_size=10 -g0.2 -N1e${SLURM_ARRAY_TASK_ID} $*
