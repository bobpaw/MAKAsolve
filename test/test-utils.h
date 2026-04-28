#ifndef _TEST_UTILS_H_
#define _TEST_UTILS_H_

#include <PCU.h>
#include <apfMDS.h>

// Create a Migration plan
apf::Migration* getParmaPlan(apf::Mesh* m, int parts);

// Load a serial mesh and partition it
apf::Mesh2* loadAndPartitionSerialMesh(char* modelFile, char* meshFile,
																			 pcu::PCU& PCU, int expectedVerts,
																			 int expectedCells);

#endif // _TEST_UTILS_H_