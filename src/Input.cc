#include <MAKAsolve/Input.h>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

namespace maka {

Input readInput(const std::string& filename) {
	Input input;

	std::ifstream file(filename);
	if (!file) {
		throw std::runtime_error("Failed to open config file: " + filename);
	}

	std::string line;

	while (std::getline(file, line)) {
		if (line.empty() || line[0] == '#') continue;

		std::istringstream iss(line);
		std::string key;
		iss >> key;

		if (key == "solver") {
			std::string value;
			std::transform(value.begin(), value.end(), value.begin(), ::tolower);
			iss >> value;
			if (value == "gpu")
				input.solver = SolverType::GPU;
			else if (value == "cpu")
				input.solver = SolverType::CPU;
			else
				throw std::runtime_error("Invalid solver: " + value);
		}

		else if (key == "kappa") {
			iss >> input.kappa;
		} else if (key == "adv_dir") {
			iss >> input.adv_dir;
		} else if (key == "adv_mag") {
			iss >> input.adv_mag;
		} else if (key == "dirichlet") {
			DirichletBC bc;
			iss >> bc.model_dim >> bc.model_tag >> bc.value;
			input.dirichletBCs.push_back(bc);
		}
	}

	return input;
}

}
