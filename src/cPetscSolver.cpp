/*
 * cPetscSolver.cpp
 *
 *  Created on: Apr 15, 2015
 *      Author: jrug001
 */

#include <iostream>
#include <unistd.h>
#include <stdlib.h>

#include "cCellMesh.h"
#include "cPetscSolver.h"

#define TEMP_SIZE 40

cPetscSolver::cPetscSolver(int argc, char **args){
    // initialize PETSc (and MPI)
	PetscInitialize(&argc, &args, PETSC_NULL, PETSC_NULL);
	PetscPushErrorHandler(PetscAbortErrorHandler, NULL);
	nonzeroguess = PETSC_FALSE;
	PetscOptionsGetBool(NULL, "-nonzero_guess", &nonzeroguess, NULL);

	// get the hostname
	char temp[TEMP_SIZE];
	gethostname(temp, TEMP_SIZE);
	host_name = temp;

	// get the MPI rank
	MPI_Comm_rank(PETSC_COMM_WORLD,&mpi_rank );
	std::cout << "<SOLVER> MPI_rank:" << mpi_rank << "  hostname:" << host_name<< std::endl;
	sync_flush();
}

cPetscSolver::~cPetscSolver(){
	VecDestroy(&x);
	VecDestroy(&b);
	MatDestroy(&A);
	KSPDestroy(&ksp);
	PetscFinalize(); // shut down PETSc (and MPI)
}

void cPetscSolver::fatal_error(std::string message){
	std::cout << message.c_str() << std::endl;
	PetscFinalize(); // shut down PETSc (and MPI)
	exit(1);
}

void cPetscSolver::init(PetscInt s, MatrixXXC &u, MatrixXXC &Amat){
	std::cout << "<SOLVER> initialising the solver..." << std::endl;
	size = s;

	// create the solution and right hand side vectors for the linear system, Ax = b
	VecCreate(PETSC_COMM_WORLD, &x);
	PetscObjectSetName((PetscObject) x, "Solution");
	VecSetSizes(x, PETSC_DECIDE, size);
	VecSetFromOptions(x);
	VecDuplicate(x, &b);

	// initialise x vector
	for(PetscInt i = 0; i < size; i++){
		PetscScalar value = u(i, 0);
		VecSetValue(x, i, value, INSERT_VALUES);
		VecAssemblyBegin(x);
		VecAssemblyEnd(x);
	}

	// create the A matrix;
	MatCreate(PETSC_COMM_WORLD, &A);
	MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, size, size);
	MatSetFromOptions(A);
	MatSetUp(A);
	MatZeroEntries(A);
	for(PetscInt i = 0; i < size; i++){
		for(PetscInt j = 0; j < size; j++){
			PetscScalar value = Amat(i, j);
			if(value != 0.0) MatSetValue(A, i, j, value, INSERT_VALUES);
		}
	}
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

	// setup linear solver
	KSPCreate(PETSC_COMM_WORLD, &ksp); // create linear solver context
	KSPSetOperators(ksp, A, A); // use the A matrix as the preconditioning matrix
	PC pc; // preconditioner context
	KSPGetPC(ksp, &pc);
	PCSetType(pc, PCJACOBI);
	//PCSetType(pc, PCILU);
	KSPSetTolerances(ksp, 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT, 250); // works for cell01m_HARMONIC_100p.msh
	//KSPSetTolerances(ksp, 1.e-9, PETSC_DEFAULT, PETSC_DEFAULT, 100); // ORIGINAL
	KSPSetFromOptions(ksp);
	KSPSetType(ksp, KSPGMRES);
    KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
}

void cPetscSolver::step(MatrixX1C &solvec, MatrixX1C &rhsvec){
	//PetscErrorCode ierr;
	for(PetscInt i = 0; i < size; i++){
		PetscScalar value = rhsvec(i);
		VecSetValue(b, i, value, INSERT_VALUES);
		VecAssemblyBegin(b);
		VecAssemblyEnd(b);
	}
	KSPSolve(ksp, b, x);
	KSPConvergedReason reason;
	KSPGetConvergedReason(ksp, &reason);
	std::cout << reason << " ";
	PetscInt its;
	KSPGetIterationNumber(ksp, &its);
	std::cout << " " << its << std::endl;
	for(PetscInt i = 0; i < size; i++){
		PetscScalar value;
		VecAssemblyBegin(x);
		VecAssemblyEnd(x);
		VecGetValues(x, 1, &i, &value);
		solvec(i) = value;
	}
}

void cPetscSolver::sync_flush(){
	PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
}
