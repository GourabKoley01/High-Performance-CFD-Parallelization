cd bld
make
cd ..
nsys profile -t nvtx,openacc,cuda,mpi -b dwarf -f true -o out mpirun -np 8 ./pravah3d.exe
#CUDA_VISIBLE_DEVICE='0' nsys profile -t nvtx,openacc,cuda -b dwarf -f true -o out mpirun -np 8 -mca coll_hcoll_enable 0 ./pravah3d.exe
