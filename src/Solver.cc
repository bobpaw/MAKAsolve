#include <MAKAsolve/Solver.h>
// #include <MAKAsolve/SparseSolver.h>
#include <cmath>
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
#include <map>
#include <mthQR.h>

namespace maka {

// extract mesh from phi, create node numbering, build BC map
Solver::Solver(apf::Field* phi, const Input& input) : phi_(phi), input_(input) {
	mesh_ = apf::getMesh(phi);

	constexpr int order = 1;
	apf::FieldShape* shape = apf::getLagrange(order);

	nbr_ = apf::numberOwnedNodes(mesh_, "nodes", shape);
	numNodes_ = apf::countNodes(nbr_);

	buildBCMap();
}

Solver::~Solver() {
	apf::destroyNumbering(nbr_);
}

// convert BCs into algebraic constraints
void Solver::buildBCMap() {
	for (auto& bc : input_.dirichletBCs) {
		// find geometric entity
		apf::ModelEntity* me = mesh_->findModelEntity(bc.model_dim, bc.model_tag);
		// get all nodes on that boundary
		apf::DynamicArray<apf::Node> nodes;
		apf::getNodesOnClosure(mesh_, me, nodes, apf::getShape(nbr_));
		// convert mesh nodes to global indices
		for (int i = 0; i < nodes.size(); i++) {
			int num = apf::getNumber(nbr_, nodes[i].entity, nodes[i].node, 0);
			dirichletMap_[num] = bc.value;
		}
	}
}

// assemble FEM system in sparse COO format
void Solver::assemble(LinearSystem& sys) {
	sys.n = numNodes_;
	sys.rhs.assign(numNodes_, 0.0);

	// temp sparse accumulator
	std::map<std::pair<int, int>, double> Kmap;
	const int dim = mesh_->getDimension();

	double kappa = input_.kappa;
	double adv_dir = input_.adv_dir;
	double adv_mag = input_.adv_mag;
	const apf::Matrix3x3 kappa_eye(kappa, 0, 0, 0, kappa, 0, 0, 0, kappa);
	apf::Vector3 adv(adv_mag * std::cos(apf::pi / 180 * adv_dir),
									 adv_mag * std::sin(apf::pi / 180 * adv_dir), 0);

	apf::MeshIterator* it = mesh_->begin(dim);
	for (apf::MeshEntity* e; (e = mesh_->iterate(it));) {
		// apf::Mesh::Type type = mesh_->getType(e);
		//  Get map from local node number to global number.
		apf::NewArray<int> ien;
		int nen = apf::getElementNumbers(nbr_, e, ien);
		apf::MeshElement* me = apf::createMeshElement(mesh_, e);
		apf::Element* el = apf::createElement(phi_, me);
		// Integration order is not equal to shape function order in all cases.
		constexpr int order = 1;
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
				if (dirichletMap_.count(ien[a])) {
					// F(ien[a]) = dirichletMap_[ien[a]];
					sys.rhs[ien[a]] = dirichletMap_[ien[a]];
					// K(ien[a], ien[a]) = 1.;
					Kmap[{ien[a], ien[a]}] = 1.0;
				} else {
					for (int b = 0; b < nen; ++b) {
						/*K(ien[a], ien[b]) +=
							-(shp_grad[a] * adv) * shp[b] * wdv
							+ (shp_grad[a] * (kappa_eye * shp_grad[b])) * wdv;*/
						Kmap[{ien[a], ien[b]}] +=
							-(shp_grad[a] * adv) * shp[b] * wdv
							+ (shp_grad[a] * (kappa_eye * shp_grad[b])) * wdv;
					}
				}
			}
		}
		apf::destroyElement(el);
		apf::destroyMeshElement(me);
	}
	mesh_->end(it);

	// convert map to COO arrays
	for (auto& entry : Kmap) {
		sys.row.push_back(entry.first.first);
		sys.col.push_back(entry.first.second);
		sys.val.push_back(entry.second);
	}
}

// solve system
void Solver::solve(LinearSystem& sys) {
	std::vector<double> solution(sys.n);

#ifdef USE_CUDA
	// gpu path
	if (input_.solver == SolverType::GPU) {
		SparseSolver solver(sys.n, sys.val.size(), sys.row.data(), sys.val.data());
		solver.solve(sys.rhs.data(), solution.data());
	} else
#endif
	{
		// cpu - dense
		mth::Matrix<double> K(sys.n, sys.n);
		mth::Vector<double> F(sys.n);

		// initialize
		for (int i = 0; i < sys.n; i++) {
			F(i) = sys.rhs[i];
			for (int j = 0; i < sys.n; j++) {
				K(i, j) = 0.0;
			}
		}
		// convert sparse to dense for cpu path
		for (size_t i = 0; i < sys.val.size(); i++) {
			K(sys.row[i], sys.col[i]) += sys.val[i];
		}
		mth::Vector<double> phi_full(sys.n);
		bool solved = mth::solveQR(K, F, phi_full);
		if (!solved) {
			throw std::runtime_error("CPU solver failed");
		}
		for (int i = 0; i < sys.n; i++) {
			solution[i] = phi_full(i);
		}
	}

	// write solution back to field
	apf::DynamicArray<apf::Node> nodes;
	apf::getNodes(nbr_, nodes);
	for (int i = 0; i < nodes.size(); i++) {
		int num = apf::getNumber(nbr_, nodes[i].entity, nodes[i].node, 0);
		apf::setScalar(phi_, nodes[i].entity, nodes[i].node, solution[num]);
	}
}

// assemble and solve
void Solver::solve() {
	LinearSystem sys;
	assemble(sys);
	solve(sys);
}

} // namespace maka
