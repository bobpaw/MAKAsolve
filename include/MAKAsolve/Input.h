#ifndef MAKASOLVE_INPUT_H
#define MAKASOLVE_INPUT_H
#include <string>
#include <vector>
#include <memory>

namespace maka {

enum class SolverType { CPU, GPU };

struct DirichletBC {
	int model_dim, model_tag;
	double value;
};

class Input {
public:
	SolverType backend_solver;

	// physics
	double kappa;
	double adv_dir;
	double adv_mag;

	// boundary conditions
	std::vector<DirichletBC> dirichletBCs;
};

using InputPtr = std::unique_ptr<Input>;

InputPtr readInput(const std::string& filename);

} // namespace maka

#endif // MAKASOLVE_INPUT_H
