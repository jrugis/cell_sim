/*
 ============================================================================
 Name        : cell_3d.cpp
 Author      : J.Rugis
 Version     :
 Copyright   : (c) 2015 J.Rugis
 Description : Cell Simulation with PETSc and MPI
 ============================================================================
 */

#define MESH_FILE "cs.msh"
#include <iostream>
#include <ctime>

#include "cCellMesh.h"
#include "cPetscSolver.h"
#include "cGeneric3dModel.h"

int main(int argc,char **args){
	cCellMesh *mesh_00;
	cPetscSolver *solver;
	cGeneric3dModel *model;
	std::clock_t start;
	double duration;

	start = std::clock(); // get the time

	// setup the solver
	//  - accepts command line arguments
	//  - note: the solver manages MPI
	solver = new cPetscSolver(argc, args);
	if(solver->mpi_rank == 0) start = std::clock();

	// read in the mesh file and display some mesh information
	mesh_00 = new cCellMesh(MESH_FILE, solver);
	mesh_00->print_info();

	// setup and run a model
	//  - reads in a model parameter file
	model = new cGeneric3dModel(mesh_00, solver);
	model->run();

	// save the results
	if(solver->mpi_rank == 0){
		model->save_results();
		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << "<MAIN> MPI_rank:0   execution time: " << duration << " sec" << std::endl;
	}

	delete model;
	delete mesh_00;
	delete solver; // this shuts down MPI

	return 0;
}
