#!/bin/bash -e
#SBATCH --job-name=bh_polaron_cont_gs # job name (shows up in the queue)
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

srun julia -tauto bosons.jl -ogscan_cont --warmup=10000 --measure=25000 -N1e7 --chunk_size=10 --load=precompute_6x6_1e7.dvec -g0.0 $*
