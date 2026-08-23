#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks=2
#SBATCH --job-name=GPU
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --partition=gpu
#SBATCH --gres=gpu:2
####SBATCH --nodelist=rmgpu013
#SBATCH --time=80:10:00
#SBATCH --mem=70GB
#SBATCH --reservation=cfd-hackathon


ulimit -s unlimited

#export PATH=/home/apps/spack/opt/spack/linux-almalinux8-skylake_avx512/gcc-8.5.0/git-2.47.0-v2svoumfx2ny7xpmh6hgpj7libh2kiwm/bin:$PATH

source /home/apps/spack/share/spack/setup-env.sh
#spack load cuda@12.4.1 /76mgrt2
spack load nvhpc@24.11 /bphfxrl
spack load git@2.47.0 /v2svoum

# only if needed export NVCOMPILER_COMM_LIBS_HOME=/home/apps/spack/opt/spack/linux-almalinux8-skylake_avx512/gcc-8.5.0/nvhpc-24.11-bphfxrlohbqt7qd5igbuvb7va2ve6znl/Linux_x86_64/24.11/comm_libs

#export LD_LIBRARY_PATH=/home/apps/spack/opt/spack/linux-almalinux8-skylake_avx512/gcc-14.2.0/gcc-runtime-14.2.0-oby3kl4y4jk2ejkret3fuq3rpyzkeydr/lib:$LD_LIBRARY_PATH

export LD_LIBRARY_PATH=/home/apps/spack/opt/spack/linux-almalinux8-skylake_avx512/gcc-14.2.0/gcc-runtime-14.2.0-oby3kl4y4jk2ejkret3fuq3rpyzkeydr/lib:$LD_LIBRARY_PATH

cd $SLURM_SUBMIT_DIR

#time ./vect_add > vector_output_gpu.txt
#runing the simutation
#echo "running the simualtion now"
#./VYOM_test.exe > log.vyomtest

#echo "Starting MPI run" 
#mpirun -np $SLURM_NTASKS ./VYOM3D_mod5.exe > log.vyomtest_mod5
#echo "Run complete" 

#nsys profile -o gpu_profile ./vect_gw

echo "profiling now"
#nsys profile -o gpu_profile_3D_vyom ./VYOM_test.exe
nsys profile \
   --trace=cuda,nvtx,mpi \
   -o gpu_profile_VYOM3D_mod5 \
   mpirun -np $SLURM_NTASKS ./VYOM3D_mod5.exe

nsys stats --report cuda_gpu_kern_sum,mpi_event_sum,nvtx_pushpop_sum gpu_profile_VYOM3D_mod5.nsys-rep

echo "profilling complete"

#echo "profiling now_txt"
##nsys stats gpu_profile_3D_vyom.nsys-rep > nsys_summary.txt
