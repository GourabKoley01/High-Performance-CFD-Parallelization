#!/bin/bash
#SBATCH --job-name=pravah_muli_g
#SBATCH --output=test_%j.output.log
#SBATCH --error=test_%j.error.log
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=60GB
#SBATCH --time=00:30:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:2
#SBATCH --account=mkvkss

### Load the modules ###
module load /home/mkvkss/rajeshr/linstall/modulefiles/nvhpc/24.5

### Run the executable ###
mpirun -np 2 ./pravah3dgpu_parallel.exe
#nsys profile -t nvtx,openacc,cuda,mpi -b dwarf -f true -o out mpirun -np 2 ./pravah3dgpu_parallel.exe

