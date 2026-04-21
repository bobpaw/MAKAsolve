#include <MAKAsolve/Input.h>
#include <MAKAsolve/Solver.h>

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

void print_exception(const std::exception& e, int level = 0);

int main(int argc, char* argv[]) {
	lion_set_verbosity(1);
	pcu::Init(&argc, &argv);
	try {
		pcu::PCU PCU;
		if (argc < 5) {
			std::cout << "USAGE: " << argv[0] << " MODEL MESH REFINEMENT [OUT.vtk]\n";
			throw 1;
		}
		char *modelFile = argv[1], *meshFile = argv[2];
		int refinement = std::stoi(argv[3]);
		char* inputFile = argv[4];
		char* vtkFile = argc > 5 ? argv[5] : NULL;

		// Initialize geometry library
		gmi_register_mesh();
		// Load mesh
		apf::Mesh2* mesh = apf::loadMdsMesh(modelFile, meshFile, &PCU);
		const int dim = mesh->getDimension();

		// Check that this is about the right file (square).
		assert(mesh->count(0) == 29);
		assert(mesh->count(2) == 40);

		// Now do refinement.
		if (refinement > 0) ma::runUniformRefinement(mesh, refinement);
		// Set finite element order.
		constexpr int order = 1;
		apf::FieldShape* shape = apf::getLagrange(order);
		apf::Field* phi = apf::createField(mesh, "phi", apf::SCALAR, shape);
		apf::zeroField(phi);

		// read input + run solver
		auto input = maka::readInput(argv[4]);
		maka::Solver solver (phi, *input);
		solver.solve();

		// Optionally plot.
		if (vtkFile) apf::writeVtkFiles(vtkFile, mesh);

		// Test for solution symmetry.
		double sum_above = 0, sum_below = 0;
		apf::MeshIterator* it = mesh->begin(0);
		for (apf::MeshEntity* vtx; (vtx = mesh->iterate(it));) {
			double value = apf::getScalar(phi, vtx, 0);
			auto x = apf::getLinearCentroid(mesh, vtx);
			if (x.y() > 0.5) sum_above += value;
			else if (x.y() < 0.5) sum_below += value;
		}
		mesh->end(it);
		double diff = std::abs(sum_above - sum_below) / mesh->count(0);
		std::cout << "above/below diff over verts: " << diff << '\n';
		if (diff > 0.01) { // Be very tolerant.
			throw std::runtime_error("failed symmetry test");
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
