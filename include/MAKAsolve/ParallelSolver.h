#pragma once

#ifdef USE_AMGX
#include <amgx_c.h>
#endif

namespace maka {

// Solves a sparse linear system Ax = b using AmgX
class ParallelSolver {
public:
	ParallelSolver(int numRows, int numNonZeros, const int* rowIndices,
								 const int* colIndices, const double* values,
								 int numDevices = 1);
	~ParallelSolver();

	void solve(const double* rhs, double* solution);

private:
	int _numRows;
	bool _distributed;

#ifdef USE_AMGX
	AMGX_config_handle _cfg;
	AMGX_resources_handle _rsrc;
	AMGX_matrix_handle _A;
	AMGX_vector_handle _b;
	AMGX_vector_handle _x;
	AMGX_solver_handle _solver;
#endif
};

} // namespace maka
