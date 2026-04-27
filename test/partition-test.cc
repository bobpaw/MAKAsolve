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

// from SCOREC core test/split.cc
apf::Migration* getPlan(apf::Mesh* m, int parts) {
	apf::Splitter* splitter = Parma_MakeRibSplitter(m);
	apf::MeshTag* weights = Parma_WeighByMemory(m);
	apf::Migration* plan = splitter->split(weights, 1.10, parts);
	apf::removeTagFromDimension(m, weights, m->getDimension());
	m->destroyTag(weights);
	delete splitter;
	return plan;
}

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
		// Load geometry model
		gmi_model* geom = gmi_load(modelFile);
		// Load unpartitioned mesh on processor 0
		apf::Mesh2* mesh = 0;
		apf::Migration* plan = 0;
		auto selfPCU = PCU.Split(PCU.Self(), PCU.Self());
		if (PCU.Self() == 0) {
			mesh = apf::loadMdsMesh(geom, meshFile, selfPCU.get());
			// Check that this is about the right file (square).
			assert(mesh->count(0) == 29);
			assert(mesh->count(2) == 40);
			plan = getPlan(mesh, PCU.Peers());
		}
		if (mesh != nullptr) mesh->switchPCU(&PCU);
		mesh = apf::repeatMdsMesh(mesh, geom, plan, PCU.Peers(), &PCU);
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
		for (apf::MeshEntity* n; (n = mesh->iterate(it));) {
			int number = apf::getNumber(gnbr, n, 0, 0);
			bool owned = mesh->isOwned(n);
			if (!owned) {
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
