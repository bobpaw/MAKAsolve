#ifndef MAKASOLVE_LINEARSOLVER_H
#define MAKASOLVE_LINEARSOLVER_H

#include <exception>
#include <memory>
#include <stdexcept>

#include <mth.h>

namespace maka {

/**
 * \brief Abstract linear solver base class.
 */
class LinearSolver {
public:
	virtual void solve(const mth::Matrix<double>& A, const mth::Vector<double>& b,
										 mth::Vector<double>& x) = 0;
	virtual ~LinearSolver() = 0;
};

class LinearException : public std::runtime_error {
public:
	using std::runtime_error::runtime_error;
};

using LinearSolverPtr = std::unique_ptr<LinearSolver>;

LinearSolverPtr makeCPULinearSolver();

} // namespace maka

#endif // MAKASOLVE_LINEARSOLVER_H
