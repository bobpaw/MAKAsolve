#pragma once

namespace maka {

// Solves a sparse linear system Ax = b on the GPU
class Solver {
 public:
  Solver(int numRows, int numNonZeros, const int* rowIndices, const int* colIndices, const double* values);
  ~Solver();
  
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