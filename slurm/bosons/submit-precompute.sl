#!/bin/bash -e
#SBATCH --job-name=bh_polaron_precompute # job name (shows up in the queue)
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
#SBATCH --array=0-1

m_array=(6 10)
N_array=(7 7)

m=${m_array[$SLURM_ARRAY_TASK_ID]}
N=${N_array[$SLURM_ARRAY_TASK_ID]}

srun julia -tauto bosons.jl -oprecompute --warmup=20000 --measure=10000 -N1e${N} -m${m} --chunk_size=10 --save=dvec_${m}x${m}_1e${N}_t0.01_u0.2.dvec -t0.01 $*
