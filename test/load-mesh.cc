#include <cassert>
#include <iostream>

#include <PCU.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_mesh.h>

int main(int argc, char* argv[]) {
	pcu::Init(&argc, &argv);
	try {
		pcu::PCU PCU;
		if (argc != 3) {
			std::cout << "USAGE: " << argv[0] << " MODEL MESH\n";
			throw 1;
		}
		char *modelFile = argv[1], *meshFile = argv[2];
		// Initialize geometry library
		gmi_register_mesh();
		// Load mesh
		apf::Mesh2* mesh = apf::loadMdsMesh(modelFile, meshFile, &PCU);
		// Check that this is about the right file (square).
		assert(mesh->count(0) == 29);
		assert(mesh->count(2) == 40);
	} catch (int r) {
		pcu::Finalize();
		return r;
	} catch (...) {
		pcu::Finalize();
		return 1;
	}
	pcu::Finalize();
	return 0;
}
