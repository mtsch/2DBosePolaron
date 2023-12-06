#!/bin/bash -e
#SBATCH --time=3-0:00:00           # Walltime (HH:MM:SS)
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

srun julia -tauto bosons.jl -otscan --warmup=20000 --measure=100000 -N1e7 --chunk_size=1000 -t${t} $*
