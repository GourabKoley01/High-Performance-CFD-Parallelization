cd bld
make
cd ..
mpirun -np 2 ./pravah3dgpu_parallel.exe
#mpirun -np 64 -mca coll_hcoll_enable 0 ./pravah3d.exe
