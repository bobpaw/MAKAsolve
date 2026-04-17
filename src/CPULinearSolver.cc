#include <MAKAsolve/LinearSolver.h>

#include <mthQR.h>

namespace maka {

class CPULinearSolver : public LinearSolver {
public:
	void solve(const mth::Matrix<double>& A, const mth::Vector<double>& b,
						 mth::Vector<double>& x) override;
	virtual ~CPULinearSolver();
};

void CPULinearSolver::solve(const mth::Matrix<double>& A, const
		mth::Vector<double>& b,
														mth::Vector<double>& x) {
	x.resize(A.cols());
	bool solved = mth::solveQR(A, b, x);
	if (!solved) throw LinearException("CPULinearSolver: failed to solve system");
}

CPULinearSolver::~CPULinearSolver() {}

LinearSolverPtr makeCPULinearSolver() {
#if __cplusplus >= 201304L
	return std::make_unique<CPULinearSolver>();
#else
	return std::unique_ptr<CPULinearSolver>(new CPULinearSolver);
#endif
}

} // namespace maka
