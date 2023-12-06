#!/bin/bash -e
#SBATCH --job-name=bh_polaron_var # job name (shows up in the queue)
#SBATCH --time=2:00:00           # Walltime (HH:MM:SS)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000                 # memory/cpu (in MB)
#SBATCH --partition=short          # specify a partition
#SBATCH --hint=nomultithread       # don't use hyperthreading
#SBATCH --error=%A.err
#SBATCH --output=%A.out
#SBATCH --mail-user=matijacufar@gmail.com
#SBATCH --mail-type=END,FAIL

srun julia -t1 variational.jl $*
