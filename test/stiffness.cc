#include <cassert>
#include <cmath>
#include <exception>
#include <iostream>
#include <string>

#include <PCU.h>
#include <apfDynamicMatrix.h>
#include <apfDynamicVector.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <lionPrint.h>
#include <gmi_mesh.h>

void print_exception(const std::exception& e, int level = 0);

int main(int argc, char* argv[]) {
	lion_set_verbosity(1);
	pcu::Init(&argc, &argv);
	try {
		pcu::PCU PCU;
		if (argc != 6) {
			std::cout << "USAGE: " << argv[0]
					<< " MODEL MESH ADV_DEG ADV_MAG KAPPA\n";
			throw 1;
		}
		char *modelFile = argv[1], *meshFile = argv[2];
		double adv_dir = std::stod(argv[3]), adv_mag = std::stod(argv[4]);
		double kappa = std::stod(argv[5]);
		const apf::Matrix3x3 kappa_eye(
			kappa, 0, 0,
			0, kappa, 0,
			0, 0, kappa
		);
		apf::Vector3 adv(
			adv_mag * std::cos(apf::pi / 180 * adv_dir),
			adv_mag * std::sin(apf::pi / 180 * adv_dir),
			0
		);
		// Initialize geometry library
		gmi_register_mesh();
		// Load mesh
		apf::Mesh2* mesh = apf::loadMdsMesh(modelFile, meshFile, &PCU);
		const int dim = mesh->getDimension();
		// Check that this is about the right file (square).
		assert(mesh->count(0) == 29);
		assert(mesh->count(2) == 40);
		// Set finite element order.
		constexpr int order = 1;
		apf::FieldShape* shape = apf::getH1Shape(order);
		apf::Field* phi = apf::createField(mesh, "phi", apf::SCALAR, shape);
		apf::zeroField(phi);
		// Count degrees of freedom.
		apf::Numbering* nbr = apf::numberOwnedNodes(mesh, "nodes", shape);
		int num_nodes = apf::countNodes(nbr);
		// Create and initialize stiffness matrix and forcing vector.
		apf::DynamicMatrix K(num_nodes, num_nodes);
		apf::DynamicVector F(num_nodes);
		for (int i = 0; i < num_nodes; ++i) {
			for (int j = 0; j < num_nodes; ++j) K(i, j) = 0.;
			F(i) = 0.;
		}
		// Set boundary conditions.
		apf::DynamicVector phi_full(num_nodes);
		phi_full.zero();
		// FIXME: actually set BCs
		// Assemble values of stiffness matrix and forcing vector.
		apf::MeshIterator* it = mesh->begin(dim);
		for (apf::MeshEntity* e; (e = mesh->iterate(it));) {
			apf::Mesh::Type type = mesh->getType(e);
			// Get map from local node number to global number.
			apf::NewArray<int> ien;
			int nen =	apf::getElementNumbers(nbr, e, ien);
			apf::MeshElement* me = apf::createMeshElement(mesh, e);
			apf::Element* el = apf::createElement(phi, me);
			// Integration order is not equal to shape function order in all cases.
			const int numIP = apf::countIntPoints(me, order);
			for (int q = 0; q < numIP; ++q) {
				apf::Vector3 xi;
				apf::getIntPoint(me, order, q, xi);
				apf::Matrix3x3 jac;
				apf::getJacobian(me, xi, jac);
				double detJ = apf::getJacobianDeterminant(jac, dim);
				double wdetJ = apf::getIntWeight(me, order, q) * detJ;
				const double phi_local = apf::getScalar(el, xi);
				apf::NewArray<double> shp;
				apf::NewArray<apf::Vector3> shp_grad;
				apf::getBF(shape, me, xi, shp);
				apf::getGradBF(shape, me, xi, shp_grad);
				for (int a = 0; a < nen; ++a) {
					F(ien[a]) += -(shp_grad[a] * adv) * phi_local * wdetJ;
					for (int b = 0; b < nen; ++b) {
						K(ien[a], ien[b]) += -(shp_grad[a] * adv) * shp[b] * wdetJ
							+ (shp_grad[a] * (kappa_eye * shp_grad[b])) * wdetJ;
					}
				}
			}
			apf::destroyElement(el);
			apf::destroyMeshElement(me);
		}
		mesh->end(it);
	} catch (int r) {
		pcu::Finalize();
		return r;
	} catch (const std::exception& e) {
		std::cerr << "ERROR: "; print_exception(e);
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
