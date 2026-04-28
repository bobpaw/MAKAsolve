#ifndef MAKASOLVE_SOLVER_H
#define MAKASOLVE_SOLVER_H

#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#include <MAKAsolve/Input.h>
#include <apfField.h>
#include <apfMesh.h>
#include <apfNumbering.h>
#include <map>
#include <utility>
#include <vector>

namespace maka {

class Solver {
public:
	Solver(apf::Field* phi, const Input& input, pcu::PCU* pcu);
	~Solver();

	// assemble and solve
	void solve();

private:
	apf::Field* phi_;
	apf::Mesh* mesh_;
	const Input& input_;

	pcu::PCU* pcu_;
	apf::GlobalNumbering* gnbr_;
	int numNodes_;
	int min_owned_;
	int max_owned_;
	int n_owned_;

	std::map<int, double> dirichletMap_;

	// convert BCs into algebraic constraints
	void buildBCMap();

	// solve in parallel with HYPRE
	void assemble(HYPRE_IJMatrix A, HYPRE_IJVector b, HYPRE_IJVector x);

	void solve(HYPRE_IJMatrix A, HYPRE_IJVector b, HYPRE_IJVector x,
						 MPI_Comm& comm);
};
} // namespace maka

#endif // MAKASOLVE_SOLVER_H
