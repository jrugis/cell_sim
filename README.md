# cell_sim
1) build your own local copy of petsc

2) git clone https://github.com/jrugis/cell_sim.git

3) build the simulation code

  a) cd cell_sim/src

  b) cp JWR_Makefile Makefile

  c) edit Makefile to point to your pets build

  d) load modules as shown in the makefile

  e) make

  f) cd ..

4) run the simulation

  a) cp JWR_run_sim.sl run_sim.sl

  b) edit run_sim.sl to point to your pets build 

  c) sbatch run_sim.sl

5) check the results

  a) diff c.bin c_REF.bin
