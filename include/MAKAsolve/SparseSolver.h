#pragma once

namespace maka {

// Solves a sparse linear system Ax = b on the GPU
class SparseSolver {
public:
	SparseSolver(int numRows, int numNonZeros, const int* rowIndices,
							 const int* colIndices, const double* values);
	~SparseSolver();

	void solve(const double* rhs, double* solution);

private:
	int _numRows;
	int _numNonZeros;

	// CSR arrays on GPU
	int* _rowOffsets;
	int* _colIndices;
	double* _values;
};

} // namespace maka
