#include <MAKAsolve/HYPRESolver.h>
#include <MAKAsolve/Input.h>
#include <MAKAsolve/Timer.h>

#include <cassert>
#include <cmath>
#include <exception>
#include <iostream>
#include <string>

#include <PCU.h>
#include <apfDynamicArray.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <gmi_mesh.h>
#include <lionPrint.h>
#include <ma.h>

#include "test-utils.h"

void print_exception(const std::exception& e, int level = 0);

int main(int argc, char* argv[]) {
	lion_set_verbosity(1);
	pcu::Init(&argc, &argv);
	try {
		pcu::PCU PCU;
		if (argc < 6) {
			std::cout << "USAGE: " << argv[0]
								<< " MODEL MESH REFINEMENT PRE_PART_REF INPUT.rc [OUT.vtk]\n";
			throw 1;
		}
		char *modelFile = argv[1], *meshFile = argv[2];
		int refinement = std::stoi(argv[3]);
		int pre_refinement = std::stoi(argv[4]); // needed for more partitions
		auto input = maka::readInput(argv[5]);
		char* vtkFile = argc > 6 ? argv[6] : NULL;

		// Timer use
		maka::Timer timer(&PCU);
		timer.prepend_info("ranks", PCU.Peers());
		timer.start_time();

		// Initialize geometry library
		gmi_register_mesh();
		// Load mesh
		apf::Mesh2* mesh = loadAndPartitionSerialMesh(modelFile, meshFile, PCU, 29,
																									40, pre_refinement);
		const int dim = mesh->getDimension();
		// Now do refinement.
		if (refinement > 0) ma::runUniformRefinement(mesh, refinement);
		// Entity counts
		int numOwnedNodes = apf::countOwned(mesh, 0);
		timer.prepend_info("num_proc_nodes", numOwnedNodes);
		int numMeshNodes = PCU.Add<double>(numOwnedNodes);
		timer.prepend_info("total_mesh_nodes", numMeshNodes);

		// Set finite element order.
		constexpr int order = 1;
		apf::Field* phi = apf::createLagrangeField(mesh, "phi", apf::SCALAR, order);
		apf::zeroField(phi);

		// read input + run solver
		timer.stop_time("setup");
		maka::HYPRESolver solver(phi, *input, &timer);
		solver.solve(&timer);

		// Optionally plot.
		if (vtkFile) apf::writeVtkFiles(vtkFile, mesh);

		// Output data
		if (PCU.Self() == 0) {
			timer.print_header_line();
			timer.print_times_line();
		}
	} catch (int r) {
		pcu::Finalize();
		return r;
	} catch (const std::exception& e) {
		std::cerr << "ERROR: ";
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
