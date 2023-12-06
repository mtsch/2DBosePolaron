#!/bin/bash -e
#SBATCH --job-name=bh_polaron_test # job name (shows up in the queue)
#SBATCH --time=2:00:00           # Walltime (HH:MM:SS)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=0                    # memory/cpu (in MB)
#SBATCH --partition=short          # specify a partition
#SBATCH --hint=nomultithread       # don't use hyperthreading
#SBATCH --error=%A
#SBATCH --output=%A
#SBATCH --mail-user=matijacufar@gmail.com
#SBATCH --mail-type=END,FAIL

srun julia -tauto bosons.jl --warmup=10000 --chunk_size=1 -g0.2 $*
