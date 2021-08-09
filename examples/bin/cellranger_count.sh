#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH --time=48:00:00   # walltime
#SBATCH --ntasks=32   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4G   # memory per CPU core
#SBATCH -J "cellranger"   # job name
#SBATCH --mail-user=dibidave@caltech.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --gid=hpc_gradinarulab
#SBATCH -o /home/davidb/slurm/slurm.%N.%j.out # STDOUT
#SBATCH -e /home/davidb/slurm/slurm.%N.%j.err # STDERR


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
