# cell_sim
1. build your own local copy of petsc
2. git clone https://github.com/jrugis/cell_sim.git
3. build the simulation code
  1. cd cell_sim/src
  2. cp JWR_Makefile Makefile
  3. edit Makefile to point to your pets build
  4. load modules as shown in the makefile
  5. make
  6. cd ..
4. run the simulation
  1. cp JWR_run_sim.sl run_sim.sl
  2. edit run_sim.sl to point to your pets build 
  3. sbatch run_sim.sl
5. check the results
  1. diff c.bin c_REF.bin
