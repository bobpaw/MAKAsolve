#pragma once
#include <map>
#include <string>
#include <vector>

namespace maka {

enum class SolverType { CPU, GPU };

struct DirichletBC {
	int model_dim;
	int model_tag;
	double value;
};

class Input {
public:
	SolverType solver;

	// physics
	double kappa;
	double adv_dir;
	double adv_mag;

	// boundary conditions
	std::vector<DirichletBC> dirichletBCs;
};

Input readInput(const std::string& filename);

} // namespace maka
