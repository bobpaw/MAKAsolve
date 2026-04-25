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

#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif

using namespace std;
static const double EPSILON = 0.0001;

void print_exception(const exception& e, int level = 0);

int main(int argc, char* argv[]) {
	lion_set_verbosity(1);
	pcu::Init(&argc, &argv);
	try {
		pcu::PCU PCU;
		if (argc < 4) {
			cout << "USAGE: " << argv[0] << " MODEL MESH INPUT [REFINEMENT]\n";
			throw 1;
		}
		char *modelFile = argv[1], *meshFile = argv[2], *inputFile = argv[3];
		int refinement = argc > 4 ? stoi(argv[4]) : 0;
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
		// Build the linear system (same for both solvers)
		maka::LinearSystem sys;
		{
			apf::zeroField(phi);
			maka::Input tmpInput = *inputPtr;
			tmpInput.backend_solver = maka::SolverType::CPU;
			maka::Solver tmpSolver(phi, tmpInput);
			tmpSolver.assemble(sys);
		}
#ifdef USE_CUDA
		// Warm up CUDA runtime before timing
		cudaFree(0);
#endif
		// Solve on CPU
		maka::Input cpuInput = *inputPtr;
		cpuInput.backend_solver = maka::SolverType::CPU;
		vector<double> cpuSol;
		apf::zeroField(phi);
		{
			maka::Solver cpuSolver(phi, cpuInput);
			auto t0 = chrono::steady_clock::now();
			cpuSolver.solve(sys);
			double cpuMs =
				chrono::duration<double, milli>(chrono::steady_clock::now() - t0)
					.count();
			cout << "CPU: " << cpuMs << " ms\n";
#ifdef USE_CUDA
			// Store CPU solution for verification against GPU later
			apf::MeshIterator* it = mesh->begin(0);
			for (apf::MeshEntity* v; (v = mesh->iterate(it));)
				cpuSol.push_back(apf::getScalar(phi, v, 0));
			mesh->end(it);
#endif
		}
#ifdef USE_CUDA
		// Solve on GPU
		maka::Input gpuInput = *inputPtr;
		gpuInput.backend_solver = maka::SolverType::GPU;
		apf::zeroField(phi);
		{
			maka::Solver gpuSolver(phi, gpuInput);
			auto t0 = chrono::steady_clock::now();
			gpuSolver.solve(sys);
			double gpuMs =
				chrono::duration<double, milli>(chrono::steady_clock::now() - t0)
					.count();
			cout << "GPU: " << gpuMs << " ms\n";
			// Compare solutions
			double diff = 0.0;
			int i = 0;
			apf::MeshIterator* it = mesh->begin(0);
			for (apf::MeshEntity* v; (v = mesh->iterate(it));)
				diff += pow(cpuSol[i++] - apf::getScalar(phi, v, 0), 2);
			mesh->end(it);
			if (sqrt(diff) < EPSILON) cout << "Results are equal\n";
			else throw runtime_error("GPU found a different solution");
		}
#else
		cout << "ERROR: USE_CUDA not defined\n";
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

void print_exception(const std::exception& e, int level) {
	std::cerr << std::string(2 * level, ' ') << e.what() << '\n';
	try {
		std::rethrow_if_nested(std::current_exception());
	} catch (const std::exception& e) {
		print_exception(e, level + 1);
	}
}
