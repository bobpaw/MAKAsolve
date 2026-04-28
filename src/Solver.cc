#include <MAKAsolve/Solver.h>
// #include <MAKAsolve/SparseSolver.h>
#include <cassert>
#include <cmath>
#include <limits.h>
#include <string>

#include <PCU.h>
#include <apfDynamicArray.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <gmi_mesh.h>
#include <lionPrint.h>
#include <ma.h>
#include <map>
#include <mthQR.h>

namespace maka {

// extract mesh from phi, create node numbering, build BC map
Solver::Solver(apf::Field* phi, const Input& input, pcu::PCU* pcu)
		: phi_(phi), input_(input), pcu_(pcu) {
	mesh_ = apf::getMesh(phi);

	constexpr int order = 1;
	apf::FieldShape* shape = apf::getLagrange(order);

	// create a global numbering
	gnbr_ = apf::makeGlobal(apf::numberOwnedNodes(mesh_, "nodes", shape));
	apf::synchronize(gnbr_);

	// find numbering information for each process
	min_owned_ = INT_MAX;
	max_owned_ = 0;
	n_owned_ = 0;
	apf::MeshIterator* it = mesh_->begin(0);
	for (apf::MeshEntity* n; (n = mesh_->iterate(it));) {
		int number = apf::getNumber(gnbr_, n, 0, 0);
		bool owned = mesh_->isOwned(n);
		if (owned) {
			min_owned_ = std::min(min_owned_, number);
			max_owned_ = std::max(max_owned_, number);
			n_owned_++;
		}
	}
	assert(max_owned_ - min_owned_ + 1 == n_owned_);
	mesh_->end(it);

	buildBCMap();
}

Solver::~Solver() {
	apf::destroyGlobalNumbering(gnbr_);
}

// convert BCs into algebraic constraints
void Solver::buildBCMap() {
	for (auto& bc : input_.dirichletBCs) {
		// find geometric entity
		apf::ModelEntity* me = mesh_->findModelEntity(bc.model_dim, bc.model_tag);
		// get all nodes on that boundary
		apf::DynamicArray<apf::Node> nodes;
		apf::getNodesOnClosure(mesh_, me, nodes, apf::getShape(gnbr_));
		// convert mesh nodes to global indices
		for (int i = 0; i < nodes.size(); i++) {
			int num = apf::getNumber(gnbr_, nodes[i].entity, nodes[i].node, 0);
			dirichletMap_[num] = bc.value;
		}
	}
}

// assemble and solve
void Solver::solve() {
	// switch between solver methods here
	if (input_.backend_solver == SolverType::CPU) {
		MPI_Comm comm;
		pcu_->DupComm(&comm);

		HYPRE_SetExecutionPolicy(HYPRE_EXEC_HOST);
		HYPRE_SetMemoryLocation(HYPRE_MEMORY_HOST);

		HYPRE_IJMatrix A;
		HYPRE_IJVector b;
		HYPRE_IJVector x;

		HYPRE_IJMatrixCreate(comm, min_owned_, max_owned_, min_owned_, max_owned_,
												 &A);
		HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
		HYPRE_IJMatrixInitialize(A);

		HYPRE_IJVectorCreate(comm, min_owned_, max_owned_, &b);
		HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
		HYPRE_IJVectorInitialize(b);

		HYPRE_IJVectorCreate(comm, min_owned_, max_owned_, &x);
		HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
		HYPRE_IJVectorInitialize(x);

		assemble(A, b, x);
		solve(A, b, x, comm);

		HYPRE_IJMatrixDestroy(A);
		HYPRE_IJVectorDestroy(b);
		HYPRE_IJVectorDestroy(x);

		MPI_Comm_free(&comm);
	}
}

// solve in parallel with HYPRE
void Solver::assemble(HYPRE_IJMatrix A, HYPRE_IJVector b, HYPRE_IJVector x) {
	double kappa = input_.kappa;
	double adv_dir = input_.adv_dir;
	double adv_mag = input_.adv_mag;
	const int dim = mesh_->getDimension();
	const apf::Matrix3x3 kappa_eye(kappa, 0, 0, 0, kappa, 0, 0, 0, kappa);
	apf::Vector3 adv(adv_mag * std::cos(apf::pi / 180 * adv_dir),
									 adv_mag * std::sin(apf::pi / 180 * adv_dir), 0);

	apf::NewArray<int> rowSizes(n_owned_);
	apf::MeshIterator* it = mesh_->begin(0);
	int r = 0;
	for (apf::MeshEntity* n; (n = mesh_->iterate(it));) {
		int number = apf::getNumber(gnbr_, n, 0, 0);
		if (mesh_->isOwned(n)) {
			apf::DynamicArray<apf::MeshEntity*> edges;
			mesh_->getAdjacent(n, 1, edges);
			rowSizes[r] = edges.getSize() + 1;
			r++;
		}
	}
	mesh_->end(it);

	HYPRE_IJMatrixSetRowSizes(A, &rowSizes[0]);

	it = mesh_->begin(dim);
	for (apf::MeshEntity* e; (e = mesh_->iterate(it));) {
		// build ien array
		apf::Downward e_nodes;
		int nen = mesh_->getDownward(e, 0, e_nodes);
		apf::NewArray<int> ien(nen);
		apf::NewArray<bool> owns(nen);
		for (int i = 0; i < nen; ++i) {
			ien[i] = apf::getNumber(gnbr_, e_nodes[i], 0, 0);
			owns[i] = mesh_->isOwned(e_nodes[i]);
		}

		// awful way to do this, sorry
		// for setting (BC nodes that this proc owns)
		int set_rows_idx = 0; // same as nrows by end of loop
		apf::NewArray<int> set_rows(nen);
		apf::NewArray<int> set_ncols(nen);
		apf::NewArray<int> set_cols(nen);
		apf::NewArray<double> set_values(nen);
		apf::NewArray<double> bc_values(nen);
		for (int i = 0; i < nen; i++) {
			set_ncols[i] = 1.0;	 // for BCs each row will only have one column
			set_values[i] = 1.0; // for BCs the one column will be 1
		}

		// for adding (non-BC nodes)
		int add_rows_idx = 0; // same as nrows by end of loop
		int add_values_idx = 0;
		apf::NewArray<int> add_rows(nen);
		apf::NewArray<int> add_ncols(nen);
		apf::NewArray<int> add_cols(nen * nen);
		apf::NewArray<double> add_values(nen * nen);
		for (int i = 0; i < nen * nen; i++) add_values[i] = 0.0;

		apf::MeshElement* me = apf::createMeshElement(mesh_, e);
		apf::Element* el = apf::createElement(phi_, me);
		// Integration order is not equal to shape function order in all cases.
		constexpr int order = 1;
		const int numIP = apf::countIntPoints(me, order);
		for (int q = 0; q < numIP; ++q) {
			apf::Vector3 xi;
			apf::getIntPoint(me, order, q, xi);
			const double dv = apf::getDV(me, xi);
			const double wdv = apf::getIntWeight(me, order, q) * dv;
			apf::NewArray<double> shp;
			apf::NewArray<apf::Vector3> shp_grad;
			apf::getShapeValues(el, xi, shp);
			apf::getShapeGrads(el, xi, shp_grad);

			add_rows_idx = 0;
			add_values_idx = 0;
			for (int a = 0; a < nen; a++) {
				bool is_bc = dirichletMap_.count(ien[a]);
				if (owns[a] && is_bc) {
					if (q == 0) {
						set_rows[set_rows_idx] = ien[a];
						set_cols[set_rows_idx] = ien[a];
						bc_values[set_rows_idx] = dirichletMap_[ien[a]];
						set_rows_idx++;
					}
				} else if (!is_bc) {
					add_rows[add_rows_idx] = ien[a];
					add_ncols[add_rows_idx] = nen;
					for (int b = 0; b < nen; b++) {
						double K_a_b_contrib =
							-(shp_grad[a] * adv) * shp[b] * wdv
							+ (shp_grad[a] * (kappa_eye * shp_grad[b])) * wdv;
						add_cols[add_values_idx] = ien[b];
						add_values[add_values_idx] += K_a_b_contrib;
						add_values_idx++;
					}
					add_rows_idx++;
				}
			}
		}
		apf::destroyElement(el);
		apf::destroyMeshElement(me);

		// indices end up being sizes after the loops
		HYPRE_IJMatrixSetValues(A, set_rows_idx, &set_ncols[0], &set_rows[0],
														&set_cols[0], &set_values[0]);
		HYPRE_IJMatrixAddToValues(A, add_rows_idx, &add_ncols[0], &add_rows[0],
															&add_cols[0], &add_values[0]);
		HYPRE_IJVectorSetValues(b, set_rows_idx, &set_rows[0], &bc_values[0]);
	}
	mesh_->end(it);

	// fill solution vector with zeros
	apf::NewArray<double> zeros(n_owned_);
	apf::NewArray<int> sol_zeros_idx(n_owned_);
	for (int i = 0; i < n_owned_; i++) {
		sol_zeros_idx[i] = min_owned_ + i;
		zeros[i] = 0.0;
	}
	HYPRE_IJVectorSetValues(x, n_owned_, &sol_zeros_idx[0], &zeros[0]);

	HYPRE_IJMatrixAssemble(A);
	HYPRE_IJVectorAssemble(b);
	HYPRE_IJVectorAssemble(x);
}

void Solver::solve(HYPRE_IJMatrix A, HYPRE_IJVector b, HYPRE_IJVector x,
									 MPI_Comm& comm) {
	// convert to inputs solver can use
	HYPRE_ParCSRMatrix parcsr_A;
	HYPRE_ParVector par_b;
	HYPRE_ParVector par_x;
	HYPRE_IJMatrixGetObject(A, (void**)&parcsr_A);
	HYPRE_IJVectorGetObject(b, (void**)&par_b);
	HYPRE_IJVectorGetObject(x, (void**)&par_x);

	// adapted from HYPRE Ex. 5
	// https://ftp.mcs.anl.gov/pub/pdetools/nightlylogs/xsdk/xsdk-configuration-tester/packages/hypre-2015.05.07/src/examples/README.html
	int num_iterations;
	double final_res_norm;
	int restart = 30;
	int modify = 1;

	/* Create solver */
	HYPRE_Solver solver, precond;
	HYPRE_ParCSRFlexGMRESCreate(comm, &solver);

	HYPRE_FlexGMRESSetKDim(solver, restart);
	HYPRE_FlexGMRESSetMaxIter(solver, 1000); /* max iterations */
	HYPRE_FlexGMRESSetTol(solver, 1e-7);		 /* conv. tolerance */
	HYPRE_FlexGMRESSetPrintLevel(solver, 2); /* print solve info */
	HYPRE_FlexGMRESSetLogging(solver, 1);		 /* needed to get run info later */

	HYPRE_BoomerAMGCreate(&precond);
	HYPRE_BoomerAMGSetPrintLevel(precond, 1);	 /* print amg solution info */
	HYPRE_BoomerAMGSetCoarsenType(precond, 8); // switch from 6 to 8
	HYPRE_BoomerAMGSetRelaxType(precond, 6);	 /* Sym G.S./Jacobi hybrid */
	HYPRE_BoomerAMGSetNumSweeps(precond, 1);
	HYPRE_BoomerAMGSetTol(precond, 0.0);	 /* conv. tolerance zero */
	HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */

	HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
														(HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup,
														precond);

	/* Now setup and solve! */
	HYPRE_ParCSRFlexGMRESSetup(solver, parcsr_A, par_b, par_x);
	HYPRE_ParCSRFlexGMRESSolve(solver, parcsr_A, par_b, par_x);

	/* Run info - needed logging turned on */
	HYPRE_FlexGMRESGetNumIterations(solver, &num_iterations);
	HYPRE_FlexGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
	if (pcu_->Self() == 0) {
		printf("\n");
		printf("Iterations = %d\n", num_iterations);
		printf("Final Relative Residual Norm = %e\n", final_res_norm);
		printf("\n");
	}

	/* Destory solver and preconditioner */
	HYPRE_ParCSRFlexGMRESDestroy(solver);
	HYPRE_BoomerAMGDestroy(precond);

	// Write solution to mesh
	apf::MeshIterator* it = mesh_->begin(0);
	apf::NewArray<double> x_local(n_owned_);
	apf::NewArray<int> x_all(n_owned_);
	for (int i = 0; i < n_owned_; i++) x_all[i] = min_owned_ + i;
	HYPRE_IJVectorGetValues(x, n_owned_, &x_all[0], &x_local[0]);
	int x_idx = 0;
	for (apf::MeshEntity* n; (n = mesh_->iterate(it));) {
		bool owned = mesh_->isOwned(n);
		if (owned) {
			apf::setScalar(phi_, n, 0, x_local[x_idx]);
			x_idx++;
		}
	}
	apf::synchronize(phi_);
	mesh_->end(it);
}

} // namespace maka