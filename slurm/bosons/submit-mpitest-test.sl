#!/bin/bash -e
#SBATCH --job-name=bh_polaron_mpitest # job name (shows up in the queue)
#SBATCH --time=0:05:00           # Walltime (HH:MM:SS)
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=0                    # memory/cpu (in MB)
#SBATCH --partition=short          # specify a partition
#SBATCH --hint=nomultithread       # don't use hyperthreading
#SBATCH --error=%A
#SBATCH --output=%A
#SBATCH --mail-user=matijacufar@gmail.com
#SBATCH --mail-type=END,FAIL

srun --mpi=pmi2 julia -tauto bosons.jl -omicro --warmup=1000 --measure=25000 --chunk_size=10 -N1e8 $*
