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
void addDirichletBCs(apf::Numbering* nbr, int model_dim, int model_tag,
										 double value, std::map<int, double>& BCs);

// Nice stable parameters that will give a vertically symmetric solution.
constexpr double adv_dir = 0, adv_mag = 1, kappa = 0.5;

int main(int argc, char* argv[]) {
	lion_set_verbosity(1);
	pcu::Init(&argc, &argv);
	try {
		pcu::PCU PCU;
		if (argc < 4) {
			std::cout << "USAGE: " << argv[0] << " MODEL MESH REFINEMENT [OUT.vtk]\n";
			throw 1;
		}
		char *modelFile = argv[1], *meshFile = argv[2];
		int refinement = std::stoi(argv[3]);
		char* vtkFile = argc > 4 ? argv[4] : NULL;
		const apf::Matrix3x3 kappa_eye(kappa, 0, 0, 0, kappa, 0, 0, 0, kappa);
		apf::Vector3 adv(adv_mag * std::cos(apf::pi / 180 * adv_dir),
										 adv_mag * std::sin(apf::pi / 180 * adv_dir), 0);
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
		// Count degrees of freedom.
		apf::Numbering* nbr = apf::numberOwnedNodes(mesh, "nodes", shape);
		int num_nodes = apf::countNodes(nbr);
		// Add Dirichlet boundary conditions.
		std::map<int, double> dirichletBCs;
		addDirichletBCs(nbr, 1, 6, 1, dirichletBCs);
		// Create and initialize stiffness matrix and forcing vector.
		mth::Matrix<double> K(num_nodes, num_nodes);
		mth::Vector<double> F(num_nodes);
		for (int i = 0; i < num_nodes; ++i) {
			for (int j = 0; j < num_nodes; ++j) K(i, j) = 0.;
			F(i) = 0.;
		}
		// Assemble values of stiffness matrix and forcing vector.
		apf::MeshIterator* it = mesh->begin(dim);
		for (apf::MeshEntity* e; (e = mesh->iterate(it));) {
			apf::Mesh::Type type = mesh->getType(e);
			// Get map from local node number to global number.
			apf::NewArray<int> ien;
			int nen = apf::getElementNumbers(nbr, e, ien);
			apf::MeshElement* me = apf::createMeshElement(mesh, e);
			apf::Element* el = apf::createElement(phi, me);
			// Integration order is not equal to shape function order in all cases.
			const int numIP = apf::countIntPoints(me, order);
			for (int q = 0; q < numIP; ++q) {
				apf::Vector3 xi;
				apf::getIntPoint(me, order, q, xi);
				const double dv = apf::getDV(me, xi);
				const double wdv = apf::getIntWeight(me, order, q) * dv;
				apf::NewArray<double> shp;
				apf::NewArray<apf::Vector3> shp_grad;
				apf::getShapeValues(el, xi, shp);
				apf::getShapeGrads(el, xi, shp_grad);
				for (int a = 0; a < nen; ++a) {
					if (dirichletBCs.count(ien[a])) {
						F(ien[a]) = dirichletBCs[ien[a]];
						K(ien[a], ien[a]) = 1.;
					} else {
						for (int b = 0; b < nen; ++b) {
							K(ien[a], ien[b]) +=
								-(shp_grad[a] * adv) * shp[b] * wdv
								+ (shp_grad[a] * (kappa_eye * shp_grad[b])) * wdv;
						}
					}
				}
			}
			apf::destroyElement(el);
			apf::destroyMeshElement(me);
		}
		mesh->end(it);
		if (refinement == 0) std::cout << "K: " << K << "F: " << F << '\n';
		// Obtain solution.
		mth::Vector<double> phi_full(num_nodes);
		bool solved = mth::solveQR(K, F, phi_full);
		if (!solved) throw std::runtime_error("failed to solve system");
		if (refinement == 0) std::cout << "phi: " << phi_full << '\n';
		apf::DynamicArray<apf::Node> nodes;
		apf::getNodes(nbr, nodes);
		for (int i = 0; i < nodes.size(); ++i) {
			int num = apf::getNumber(nbr, nodes[i].entity, nodes[i].node, 0);
			apf::setScalar(phi, nodes[i].entity, nodes[i].node, phi_full(i));
		}
		apf::destroyNumbering(nbr);
		// Optionally plot.
		if (vtkFile) apf::writeVtkFiles(vtkFile, mesh);
		// Test for solution symmetry.
		double sum_above = 0, sum_below = 0;
		it = mesh->begin(0);
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

void addDirichletBCs(apf::Numbering* nbr, int model_dim, int model_tag,
										 double value, std::map<int, double>& BCs) {
	apf::Mesh* m = apf::getMesh(nbr);
	apf::ModelEntity* me = m->findModelEntity(model_dim, model_tag);
	apf::DynamicArray<apf::Node> bdry_nodes;
	apf::getNodesOnClosure(m, me, bdry_nodes, apf::getShape(nbr));
	for (int i = 0; i < bdry_nodes.size(); ++i) {
		int num = apf::getNumber(nbr, bdry_nodes[i].entity, bdry_nodes[i].node, 0);
		BCs[num] = value;
	}
}
