/*
 * cGeneric3dModel.h
 *
 *  Created on: 23/10/2015
 *      Author: jrugis
 */

#ifndef CGENERIC3DMODEL_H_
#define CGENERIC3DMODEL_H_

#include <Eigen/Dense>
#include "cCellMesh.h"
#include "cPetscSolver.h"

#define REF_MASS_SIZE 4   // reference mass dimension
#define VARIABLES 3 // the number of model variables

// the generic 3D model parameters
#define MODEL_FILE "cs.dat"
enum model_parameters{delt, tend, \
	IPRdn, IPRdf, IPRmin, PLCds, PLCdl, \
	c0, ip0, h0, ct, gama, \
	kIPR, Kc, Kp, Ki, tau, Dc, \
	VS, KS, kleak, \
	VPLC, Vdeg, K3K, Dp, \
	PCOUNT};

enum model_node_values{IPR_n, PLC_n, MODELNCOUNT};            // node spatial factors
enum model_element_values{VOL_e, IPR_e, PLC_e, MODELECOUNT};  // element volume and spatial factors

// some convenience typedefs
typedef Eigen::Matrix<tCalcs, Eigen::Dynamic, Eigen::Dynamic> MatrixXXC;
typedef Eigen::Matrix<tCalcs, Eigen::Dynamic, 1> MatrixX1C;
typedef Eigen::Array<tCalcs, Eigen::Dynamic, 1> ArrayX1C;
typedef Eigen::Array<tCalcs, 1, VARIABLES> Array1VC;
typedef Eigen::Array<tCalcs, REF_MASS_SIZE, REF_MASS_SIZE> ArrayRefMass;

class cGeneric3dModel {
public:
	cGeneric3dModel(cCellMesh *mesh, cPetscSolver *solver);
	virtual ~cGeneric3dModel();
	void run();
	void save_results();

	MatrixXXC Amat, mass, u; // A, mass and solution matrices

private:
	void get_parameters();
	void init_u();
	MatrixX1C make_load(long i);
	ArrayRefMass make_ref_mass();
	Array1VC getbodyreactions(tCalcs c, tCalcs ip, tCalcs h, tCalcs ipr_f, tCalcs plc_f);
	tCalcs getboundaryflux(tCalcs c);
	void make_matrices();
	void load_node_data(std::string file_name, int dindex);
	void save_matrix(std::string file_name, MatrixXXC mat);

	cCellMesh *mesh;
	cPetscSolver *solver;
	tCalcs p[PCOUNT]; // the model parameters array
	long numt; // number of time steps
	Eigen::Array<tCalcs, Eigen::Dynamic, MODELNCOUNT> node_data;
	Eigen::Array<tCalcs, Eigen::Dynamic, MODELECOUNT> element_data;
};

#endif /* CGENERIC3DMODEL_H_ */
