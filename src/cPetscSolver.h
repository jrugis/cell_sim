/*
 * cPetscSolver.h
 *
 *  Created on: Apr 15, 2015
 *      Author: jrug001
 */

#ifndef CPETSCSOLVER_H_
#define CPETSCSOLVER_H_

typedef double tCalcs;

#include <string>
#include <Eigen/Dense>
#include <petscksp.h>

typedef Eigen::Matrix<tCalcs, Eigen::Dynamic, Eigen::Dynamic> MatrixXXC;
typedef Eigen::Matrix<tCalcs, Eigen::Dynamic, 1> MatrixX1C;

class cPetscSolver {
public:
	cPetscSolver(int argc, char **args);
	virtual ~cPetscSolver();
	void step(MatrixX1C &solvec, MatrixX1C &rhsvec);
	void init(PetscInt s, MatrixXXC &u, MatrixXXC &Amat);
	void sync_flush();
	void fatal_error(std::string message);

	int mpi_rank;
	std::string host_name;

private:
	Vec x, b; // approximate solution, RHS
	Mat A; // linear system matrix
	KSP ksp; // linear solver context
	PetscInt size;
	PetscErrorCode ierr;
	PetscBool nonzeroguess;
};

#endif /* CPETSCSOLVER_H_ */
