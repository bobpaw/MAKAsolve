#include <MAKAsolve/ParallelSolver.h>
#include <mpi.h>
#include <stdexcept>
#include <vector>

namespace maka {

// AmgX solver configuration
static const char* AMGX_CONFIG = "config_version=2, "
                                 "solver(main)=PBICGSTAB, "
                                 "main:preconditioner(amg)=AMG, "
                                 "amg:algorithm=CLASSICAL, "
                                 "amg:max_iters=1, "
                                 "amg:cycle=V, "
                                 "main:max_iters=100, "
                                 "main:convergence=RELATIVE_INI_CORE, ";

// Convert COO from Assemble() to CSR (required for AmgX)
static std::vector<int> cooToRowOffsets(int numRows, int numNonZeros,
																				const int* rowIndices) {
	std::vector<int> rowOffsets(numRows + 1, 0);
	for (int i = 0; i < numNonZeros; i++) rowOffsets[rowIndices[i] + 1]++;
	for (int i = 0; i < numRows; i++) rowOffsets[i + 1] += rowOffsets[i];
	return rowOffsets;
}

// Get local CSR data for a range of rows
static void getLocalCSR(int rowStart, int rowEnd, int numNonZeros,
												const int* rowIndices, const int* colIndices,
												const double* values, std::vector<int>& localRowOffsets,
												std::vector<int>& localColIndices,
												std::vector<double>& localValues) {
	int localRows = rowEnd - rowStart;
	localRowOffsets.assign(localRows + 1, 0);
	// Count nonzeros in each local row
	for (int i = 0; i < numNonZeros; i++) {
		int r = rowIndices[i];
		if (r >= rowStart && r < rowEnd) localRowOffsets[r - rowStart + 1]++;
	}
	// Prefix sum to get row start positions
	for (int i = 0; i < localRows; i++)
		localRowOffsets[i + 1] += localRowOffsets[i];
	// Fill local col/val arrays using a per-row position counter
	int localNnz = localRowOffsets[localRows];
	localColIndices.resize(localNnz);
	localValues.resize(localNnz);
	std::vector<int> pos(localRows, 0);
	for (int i = 0; i < numNonZeros; i++) {
		int r = rowIndices[i];
		if (r >= rowStart && r < rowEnd) {
			int lr = r - rowStart;
			int idx = localRowOffsets[lr] + pos[lr]++;
			localColIndices[idx] = colIndices[i];
			localValues[idx] = values[i];
		}
	}
}

ParallelSolver::ParallelSolver(int numRows, int numNonZeros,
															 const int* rowIndices, const int* colIndices,
															 const double* values, int numDevices)
		: _numRows(numRows) {
	int rank = 0, size = 1;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	_distributed = (size > 1);

	AMGX_initialize();
	AMGX_config_create(&_cfg, AMGX_CONFIG);

	if (_distributed) {
		// Split available GPUs among all processes
		int device = rank % numDevices;
		MPI_Comm comm = MPI_COMM_WORLD;
		AMGX_resources_create(&_rsrc, _cfg, &comm, 1, &device);
	} else {
		// Single process, assign all available GPUs.
		std::vector<int> devices(numDevices);
		for (int i = 0; i < numDevices; i++) devices[i] = i;
		AMGX_resources_create(&_rsrc, _cfg, NULL, numDevices, devices.data());
	}

	// Create matrix and vectors in AmgX format
	AMGX_matrix_create(&_A, _rsrc, AMGX_mode_dDDI);
	AMGX_vector_create(&_b, _rsrc, AMGX_mode_dDDI);
	AMGX_vector_create(&_x, _rsrc, AMGX_mode_dDDI);
	AMGX_solver_create(&_solver, _rsrc, AMGX_mode_dDDI, _cfg);

	if (_distributed) {
		// Each rank owns a block of rows.
		int rowStart = rank * numRows / size;
		int rowEnd = (rank + 1) * numRows / size;
		int localRows = rowEnd - rowStart;

		std::vector<int> localRowOffsets, localColIndices;
		std::vector<double> localValues;
		getLocalCSR(rowStart, rowEnd, numNonZeros, rowIndices, colIndices, values,
								localRowOffsets, localColIndices, localValues);

		// Build partition vector (Used this to route halo communication)
		std::vector<int> partition(numRows);
		for (int i = 0; i < numRows; i++) partition[i] = i * size / numRows;

		// Build distribution handle with partition data.
		AMGX_distribution_handle dist;
		AMGX_distribution_create(&dist, _cfg);
		AMGX_distribution_set_partition_data(dist, AMGX_DIST_PARTITION_VECTOR,
																				 partition.data());
		AMGX_distribution_set_32bit_colindices(dist, 1);

		// Each rank uploads only its local rows
		// and AmgX handles the distribution internally
		if (AMGX_matrix_upload_distributed(
					_A, numRows, localRows, (int)localValues.size(), 1, 1,
					localRowOffsets.data(), localColIndices.data(), localValues.data(),
					NULL, dist)
				!= AMGX_RC_OK)
			throw std::runtime_error("AmgX distributed matrix upload failed");

		AMGX_distribution_destroy(dist);
	} else {
		// Single process
		std::vector<int> rowOffsets =
			cooToRowOffsets(numRows, numNonZeros, rowIndices);
		if (AMGX_matrix_upload_all(_A, numRows, numNonZeros, 1, 1,
															 rowOffsets.data(), colIndices, values, NULL)
				!= AMGX_RC_OK)
			throw std::runtime_error("AmgX matrix upload failed");
	}

	// Setup solver with matrix
	AMGX_solver_setup(_solver, _A);
}

ParallelSolver::~ParallelSolver() {
	AMGX_solver_destroy(_solver);
	AMGX_vector_destroy(_x);
	AMGX_vector_destroy(_b);
	AMGX_matrix_destroy(_A);
	AMGX_resources_destroy(_rsrc);
	AMGX_config_destroy(_cfg);
	AMGX_finalize();
}

void ParallelSolver::solve(const double* rhs, double* solution) {
	int rank = 0, size = 1;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (_distributed) {
		int rowStart = rank * _numRows / size;
		int rowEnd = (rank + 1) * _numRows / size;
		int localRows = rowEnd - rowStart;

		// Each rank slices its portion of the global rhs.
		std::vector<double> localRhs(rhs + rowStart, rhs + rowEnd);
		std::vector<double> localX(localRows, 0.0);

		AMGX_vector_upload(_b, localRows, 1, localRhs.data());
		AMGX_vector_upload(_x, localRows, 1, localX.data());

		// Handles the distributed solve internally
		AMGX_solver_solve(_solver, _b, _x);

		AMGX_SOLVE_STATUS status;
		AMGX_solver_get_status(_solver, &status);
		if (status != AMGX_SOLVE_SUCCESS)
			throw std::runtime_error("AmgX solver did not converge");

		// Download local solution and gather to all ranks.
		std::vector<double> localSol(localRows);
		AMGX_vector_download(_x, localSol.data());

		// Reassemble the full solution on all ranks
		std::vector<int> recvCounts(size), displs(size);
		for (int r = 0; r < size; r++) {
			recvCounts[r] = (r + 1) * _numRows / size - r * _numRows / size;
			displs[r] = r * _numRows / size;
		}
		MPI_Allgatherv(localSol.data(), localRows, MPI_DOUBLE, solution,
									 recvCounts.data(), displs.data(), MPI_DOUBLE,
									 MPI_COMM_WORLD);
	} else {
		// Single process
		std::vector<double> x0(_numRows, 0.0);
		AMGX_vector_upload(_b, _numRows, 1, rhs);
		AMGX_vector_upload(_x, _numRows, 1, x0.data());
		AMGX_solver_solve(_solver, _b, _x);

		// Solve directly
		AMGX_SOLVE_STATUS status;
		AMGX_solver_get_status(_solver, &status);
		if (status != AMGX_SOLVE_SUCCESS)
			throw std::runtime_error("AmgX solver did not converge");

		AMGX_vector_download(_x, solution);
	}
}

} // namespace maka
