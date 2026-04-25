#include <MAKAsolve/Input.h>
#include <cassert>
#include <iostream>

int main(int argc, char* argv[]) {
	if (argc < 2) {
		std::cerr << "USAGE: " << argv[0] << " <input_file>\n";
		return 1;
	}

	maka::InputPtr input = maka::readInput(argv[1]);

	// physics
	assert(input->kappa == 0.5);
	assert(input->adv_dir == 0.0);
	assert(input->adv_mag == 1.0);

	// solver
	assert(input->backend_solver == maka::SolverType::GPU);

	// BCs
	assert(input->dirichletBCs.size() == 4);
	assert(input->dirichletBCs[0].model_dim == 1);
	assert(input->dirichletBCs[0].model_tag == 6);
	assert(input->dirichletBCs[0].value == 1.0);
	assert(input->dirichletBCs[1].model_tag == 2);
	assert(input->dirichletBCs[2].model_tag == 4);
	assert(input->dirichletBCs[3].model_tag == 8);

	return 0;
}
