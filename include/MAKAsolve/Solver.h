#pragma once

#include <MAKAsolve/Input.h>
#include <apfField.h>
#include <apfMesh.h>
#include <apfNumbering.h>
#include <map>
#include <utility>
#include <vector>

namespace maka {

struct LinearSystem {
	// number of rows
	int n;
	// COO sparse matrix storage
	std::vector<int> row, col;
	std::vector<double> val, rhs;
};

class Solver {
public:
	Solver(apf::Field* phi, const Input& input);
	~Solver();

	// assemble, output goes into sparse linear system
	void assemble(LinearSystem& sys);
	// solve a pre assembled system and write back to phi field
	void solve(LinearSystem& sys);
	// assemble and solve
	void solve();

private:
	apf::Field* phi_;
	apf::Mesh* mesh_;
	const Input& input_;

	apf::Numbering* nbr_;
	int numNodes_;

	std::map<int, double> dirichletMap_;

	// convert BCs into algebraic constraints
	void buildBCMap();
};

} // namespace maka
