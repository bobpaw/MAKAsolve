#include <MAKAsolve/Input.h>
#include <MAKAsolve/Solver.h>
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
#include <apfNumbering.h>
#include <apfShape.h>
#include <gmi_mesh.h>
#include <lionPrint.h>
#include <ma.h>
#include <mthQR.h>

#include "test_utils.h"

void print_exception(const std::exception& e, int level = 0);

int main(int argc, char* argv[]) {
	lion_set_verbosity(1);
	pcu::Init(&argc, &argv);
	try {
		pcu::PCU PCU;
		if (argc < 5) {
			std::cout << "USAGE: " << argv[0]
								<< " MODEL MESH REFINEMENT INPUT.rc [OUT.vtk]\n";
			throw 1;
		}
		char *modelFile = argv[1], *meshFile = argv[2];
		int refinement = std::stoi(argv[3]);
		char* inputFile = argv[4];
		char* vtkFile = argc > 5 ? argv[5] : NULL;

		// Timer use
		maka::Timer timer;
		timer.prepend_info("ranks", PCU.Peers());
		timer.start_time(&PCU);

		// Initialize geometry library
		gmi_register_mesh();
		// Load mesh
		apf::Mesh2* mesh =
			loadAndPartitionSerialMesh(modelFile, meshFile, PCU, 29, 40);
		const int dim = mesh->getDimension();
		// Now do refinement.
		if (refinement > 0) ma::runUniformRefinement(mesh, refinement);
		// Set finite element order.
		constexpr int order = 1;
		apf::Field* phi = apf::createLagrangeField(mesh, "phi", apf::SCALAR, order);
		apf::zeroField(phi);

		// read input + run solver
		auto input = maka::readInput(argv[4]);
		timer.stop_time("setup", &PCU);
		maka::Solver solver(phi, *input, &PCU, &timer);
		solver.solve(&timer);

		// Optionally plot.
		if (vtkFile) apf::writeVtkFiles(vtkFile, mesh);

		// Test for solution symmetry.
		double sum_above = 0, sum_below = 0;
		apf::MeshIterator* it = mesh->begin(0);
		for (apf::MeshEntity* vtx; (vtx = mesh->iterate(it));) {
			if (!mesh->isOwned(vtx)) continue;
			double value = apf::getScalar(phi, vtx, 0);
			auto x = apf::getLinearCentroid(mesh, vtx);
			if (x.y() > 0.5) sum_above += value;
			else if (x.y() < 0.5) sum_below += value;
		}
		mesh->end(it);
		sum_above = PCU.Add<double>(sum_above);
		sum_below = PCU.Add<double>(sum_below);
		double diff = std::abs(sum_above - sum_below) / mesh->count(0);
		std::cout << "above/below diff over verts: " << diff << '\n';
		if (diff > 0.01) { // Be very tolerant.
			throw std::runtime_error("failed symmetry test");
		}

		// Example output data (printing in this format will make processing easier)
		if (PCU.Self() == 0) {
			if (PCU.Peers() == 1) timer.print_header_line();
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
