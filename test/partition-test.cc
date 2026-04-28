#include <cassert>
#include <iostream>

#include <PCU.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <gmi_mesh.h>
#include <ma.h>
#include <parma.h>

#include "test_utils.h"

int main(int argc, char* argv[]) {
	pcu::Init(&argc, &argv);
	try {
		pcu::PCU PCU;
		if (argc < 4 && PCU.Self() == 0) {
			std::cout << "USAGE: " << argv[0] << " MODEL MESH REFINEMENT [OUT.vtk]\n";
			throw 1;
		}
		char *modelFile = argv[1], *meshFile = argv[2];
		int refinement = std::stoi(argv[3]);
		char* vtkFile = argc > 4 ? argv[4] : NULL;
		// Initialize geometry library
		gmi_register_mesh();
		// Load mesh
		apf::Mesh2* mesh =
			loadAndPartitionSerialMesh(modelFile, meshFile, PCU, 29, 40);
		assert(mesh && "mesh part does not exist on processor");

		if (refinement > 0) ma::runUniformRefinement(mesh, refinement);

		// check numberings (make sure I understand them correctly)
		constexpr int order = 1;
		apf::FieldShape* shape = apf::getLagrange(order);
		apf::GlobalNumbering* gnbr =
			apf::makeGlobal(apf::numberOwnedNodes(mesh, "nodes", shape));
		apf::synchronize(gnbr);

		apf::MeshIterator* it = mesh->begin(0);
		int min_owned = ~(1 << 31);
		int max_owned = 0;
		for (apf::MeshEntity* n; (n = mesh->iterate(it));) {
			int number = apf::getNumber(gnbr, n, 0, 0);
			bool owned = mesh->isOwned(n);
			if (owned) {
				min_owned = std::min(min_owned, number);
				max_owned = std::max(max_owned, number);
			}
		}
		mesh->end(it);
		it = mesh->begin(0);
		int num_not_owned = 0;
		int last_owned_id = -1;
		for (apf::MeshEntity* n; (n = mesh->iterate(it));) {
			int number = apf::getNumber(gnbr, n, 0, 0);
			bool owned = mesh->isOwned(n);
			if (owned) {
				if (last_owned_id < 0) last_owned_id = number - 1;
				// if (PCU.Self() == 1) printf("%d vs last number %d\n", number,
				// last_owned_id);
				assert(number - last_owned_id == 1);
				last_owned_id = number;
			} else {
				num_not_owned++;
				assert(number < min_owned || number > max_owned);
			}
		}
		printf("proc %d min owned id %d max owned id %d num not owned %d\n",
					 PCU.Self(), min_owned, max_owned, num_not_owned);
		mesh->end(it);

		apf::destroyGlobalNumbering(gnbr);

		if (vtkFile) apf::writeVtkFiles(vtkFile, mesh);

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
