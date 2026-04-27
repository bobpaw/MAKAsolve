#include <MAKAsolve/Input.h>
#include <MAKAsolve/Solver.h>
#include <chrono>
#include <cmath>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include <PCU.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <gmi_mesh.h>
#include <lionPrint.h>
#include <ma.h>

#ifdef USE_AMGX
#include <cuda_runtime.h>
#endif

using namespace std;
static const bool DEBUG_COMPARE = true; // To verify correct solution
static const double EPSILON = 0.0001;

void print_exception(const exception& e, int level = 0);

int main(int argc, char* argv[]) {
	lion_set_verbosity(1);
	pcu::Init(&argc, &argv);
	try {
		pcu::PCU PCU;
		if (argc < 6) {
			cout << "USAGE: " << argv[0]
					 << " MODEL MESH INPUT NUM_NODES NUM_GPUS_PER_NODE [REFINEMENT]\n";
			throw 1;
		}
		char *modelFile = argv[1], *meshFile = argv[2], *inputFile = argv[3];
		int numNodes = stoi(argv[4]), numGPUsPerNode = stoi(argv[5]),
				refinement = argc > 6 ? stoi(argv[6]) : 0;
		// Initialize geometry library and load mesh
		gmi_register_mesh();
		apf::Mesh2* mesh = apf::loadMdsMesh(modelFile, meshFile, &PCU);
		if (refinement > 0) ma::runUniformRefinement(mesh, refinement);
		// Read physics and BC parameters
		maka::InputPtr inputPtr = maka::readInput(inputFile);
		// Create phi field
		constexpr int order = 1;
		apf::FieldShape* shape = apf::getLagrange(order);
		apf::Field* phi = apf::createField(mesh, "phi", apf::SCALAR, shape);
		// Assemble the linear system
		maka::LinearSystem sys;
		{
			apf::zeroField(phi);
			maka::Input tmpInput = *inputPtr;
			tmpInput.backend_solver = maka::SolverType::CPU;
			maka::Solver tmpSolver(phi, tmpInput);
			tmpSolver.assemble(sys);
		}
		int rank = 0;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if (rank == 0) {
			cout << "nodes=" << numNodes << " gpus_per_node=" << numGPUsPerNode
					 << " total_gpus=" << numNodes * numGPUsPerNode << "\n\n";
		}
#ifdef USE_AMGX
		// Warm up CUDA runtime before timing
		cudaFree(0);
#endif
		// Solve on CPU for verification if debug mode is on
		vector<double> cpuSol;
		if (DEBUG_COMPARE) {
			maka::Input cpuInput = *inputPtr;
			cpuInput.backend_solver = maka::SolverType::CPU;
			apf::zeroField(phi);
			{
				maka::Solver cpuSolver(phi, cpuInput);
				auto t0 = chrono::steady_clock::now();
				cpuSolver.solve(sys);
				double cpuMs =
				chrono::duration<double, milli>(chrono::steady_clock::now() - t0)
					.count();
				cout << "CPU: " << cpuMs << " ms\n";
			}
			apf::MeshIterator* it = mesh->begin(0);
			for (apf::MeshEntity* v; (v = mesh->iterate(it));)
				cpuSol.push_back(apf::getScalar(phi, v, 0));
			mesh->end(it);
		}
#ifdef USE_AMGX
		// Solve with AmgX
		maka::Input amgxInput = *inputPtr;
		amgxInput.backend_solver = maka::SolverType::AMGX;
		amgxInput.numDevices = numGPUsPerNode;
		apf::zeroField(phi);
		{
			maka::Solver amgxSolver(phi, amgxInput);
			// Barrier to synchronize before timing
			MPI_Barrier(MPI_COMM_WORLD);
			auto t0 = chrono::steady_clock::now();
			amgxSolver.solve(sys);
			MPI_Barrier(MPI_COMM_WORLD);
			if (rank == 0) {
				double ms =
					chrono::duration<double, milli>(chrono::steady_clock::now() - t0)
						.count();
				cout << "AmgX: " << ms << " ms\n";
			}
		}
		// Compare against CPU solution if debug mode is on
		if (DEBUG_COMPARE && rank == 0) {
			double diff = 0.0;
			int i = 0;
			apf::MeshIterator* it = mesh->begin(0);
			for (apf::MeshEntity* v; (v = mesh->iterate(it));)
				diff += pow(cpuSol[i++] - apf::getScalar(phi, v, 0), 2);
			mesh->end(it);
			if (sqrt(diff) < EPSILON) cout << "Results are equal\n";
			else throw runtime_error("AmgX found a different solution");
		}
#else
		cout << "ERROR: USE_AMGX not defined\n";
#endif
	} catch (int r) {
		pcu::Finalize();
		return r;
	} catch (const exception& e) {
		cerr << "ERROR: ";
		print_exception(e);
		pcu::Finalize();
		return 1;
	} catch (...) {
		pcu::Finalize();
		return 1;
	}
	pcu::Finalize();
	return 0;
}

void print_exception(const exception& e, int level) {
	cerr << string(2 * level, ' ') << e.what() << '\n';
	try {
		rethrow_if_nested(current_exception());
	} catch (const exception& nested) {
		print_exception(nested, level + 1);
	}
}
