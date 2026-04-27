#include "test_utils.h"

#include <cassert>
#include <gmi_mesh.h>
#include <parma.h>

// Create a Migration plan
apf::Migration* getParmaPlan(apf::Mesh* m, int parts) {
	apf::Splitter* splitter = Parma_MakeRibSplitter(m);
	apf::MeshTag* weights = Parma_WeighByMemory(m);
	apf::Migration* plan = splitter->split(weights, 1.10, parts);
	apf::removeTagFromDimension(m, weights, m->getDimension());
	m->destroyTag(weights);
	delete splitter;
	return plan;
}

// Load a serial mesh and partition it
apf::Mesh2* loadAndPartitionSerialMesh(char* modelFile, char* meshFile,
																			 pcu::PCU& PCU, int expectedVerts,
																			 int expectedCells) {
	gmi_register_mesh();
	// Load geometry model
	gmi_model* geom = gmi_load(modelFile);
	// Load unpartitioned mesh on processor 0
	apf::Mesh2* mesh = 0;
	apf::Migration* plan = 0;
	auto selfPCU = PCU.Split(PCU.Self(), PCU.Self());
	if (PCU.Self() == 0) {
		mesh = apf::loadMdsMesh(geom, meshFile, selfPCU.get());
		assert(mesh->count(0) == expectedVerts);
		assert(mesh->count(2) == expectedCells);
		plan = getParmaPlan(mesh, PCU.Peers());
	}
	if (mesh != nullptr) mesh->switchPCU(&PCU);
	return apf::repeatMdsMesh(mesh, geom, plan, PCU.Peers(), &PCU);
}