#!/bin/sh

#SBATCH -J jobname
#SBATCH -p BP10
#SBATCH -N 1
#SBATCH -n 1

sample=list

sh run2_sentieon.sh $sample $SLURM_NPROCS


