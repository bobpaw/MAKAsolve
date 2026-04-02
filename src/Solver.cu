#include <MAKAsolve/Solver.h>
#include <cusparse.h>
#include <cusolverSp.h>
#include <cuda_runtime.h>

namespace maka {

// Convert COO to CSR and store on GPU
// CSR format is necessary for cusolverSpDcsrlsvqr
Solver::Solver(int numRows, int numNonZeros, const int* rowIndices, const int* colIndices, const double* values)
      : _numRows(numRows), _numNonZeros(numNonZeros) {
  
  // Allocate GPU memory
  cudaMalloc(&_rowOffsets, (_numRows + 1) * sizeof(int));
  cudaMalloc(&_colIndices, _numNonZeros * sizeof(int));
  cudaMalloc(&_values, _numNonZeros * sizeof(double));

  // Copy values to GPU for CSR format
  cudaMemcpy(_colIndices, colIndices, _numNonZeros * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(_values, values, _numNonZeros * sizeof(double), cudaMemcpyHostToDevice);

  // Temporarily copy row indices to GPU
  int* rowIndicesGpu;
  cudaMalloc(&rowIndicesGpu, _numNonZeros * sizeof(int));
  cudaMemcpy(rowIndicesGpu, rowIndices, _numNonZeros * sizeof(int), cudaMemcpyHostToDevice);

  // Convert COO row indices to CSR row offsets
  cusparseHandle_t sparseHandle;
  cusparseCreate(&sparseHandle);
  cusparseXcoo2csr(sparseHandle, rowIndicesGpu, _numNonZeros, _numRows,
                   _rowOffsets, CUSPARSE_INDEX_BASE_ZERO);
  
  cusparseDestroy(sparseHandle);
  cudaFree(rowIndicesGpu);
}

Solver::~Solver() {
  cudaFree(_rowOffsets);
  cudaFree(_colIndices);
  cudaFree(_values);
}

// Solve Ax = b using QR factorization
void Solver::solve(const double* rhs, double* solution) {
  cusolverSpHandle_t solverHandle;
  cusolverSpCreate(&solverHandle);

  cusparseMatDescr_t matDescr;
  cusparseCreateMatDescr(&matDescr);
  cusparseSetMatType(matDescr, CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(matDescr, CUSPARSE_INDEX_BASE_ZERO);

  // Allocate GPU memory for rhs and solution
  double* rhsGpu;
  double* solutionGpu;
  cudaMalloc(&rhsGpu, _numRows * sizeof(double));
  cudaMalloc(&solutionGpu, _numRows * sizeof(double));
  cudaMemcpy(rhsGpu, rhs, _numRows * sizeof(double), cudaMemcpyHostToDevice);

  // Solve Ax = b using QR factorization
  int singularRow;
  cusolverSpDcsrlsvqr(solverHandle, _numRows, _numNonZeros, matDescr,
                      _values, _rowOffsets, _colIndices,
                      rhsGpu, 1e-12, 0, solutionGpu, &singularRow);
  
  // Matrix is singular
  if (singularRow >= 0) {
    cudaFree(rhsGpu);
    cudaFree(solutionGpu);
    cusparseDestroyMatDescr(matDescr);
    cusolverSpDestroy(solverHandle);

    return; // Return empty solution
  }

  // Copy solution back to CPU
  cudaMemcpy(solution, solutionGpu, _numRows * sizeof(double), cudaMemcpyDeviceToHost);
  
  cudaFree(rhsGpu);
  cudaFree(solutionGpu);
  cusparseDestroyMatDescr(matDescr);
  cusolverSpDestroy(solverHandle);
}

} // namespace maka