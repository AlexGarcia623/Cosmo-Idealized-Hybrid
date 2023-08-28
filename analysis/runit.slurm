#!/bin/bash					
#SBATCH --job-name=mov
#SBATCH --mem-per-cpu=5000mb
#SBATCH --time=0-3:00:00
#SBATCH --mail-user=j.rose@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=hpg2-compute
#SBATCH --ntasks=10
##SBATCH --dependency=singleton
#SBATCH --qos=paul.torrey
##SBATCH --account=astronomy-dept
##SBATCH --qos=astronomy-dept
#SBATCH --output=plots/logs/output_%j.out
#SBATCH --error=plots/logs/error_%j.err

module purge
module load intel/2018.1.163 gsl/2.4 openmpi/3.1.2 python3
module list

mpirun -n 10 python make_movie.py 
