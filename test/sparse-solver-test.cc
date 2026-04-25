#include <MAKAsolve/SparseSolver.h>
#include <cassert>
#include <cmath>
#include <iostream>

using namespace std;

static const double EPSILON = 0.0001;

static void assertSolution(const double* solution, const double* expected,
													 int n) {
	for (int i = 0; i < n; i++)
		assert(fabs(solution[i] - expected[i]) <= EPSILON);
}

static void runTest(int numRows, int numNonZeros, int* rowIndices,
										int* colIndices, double* values, double* rhs,
										double* expected) {

	double* solution = new double[numRows]();

	maka::SparseSolver solver(numRows, numNonZeros, rowIndices, colIndices,
														values);
	solver.solve(rhs, solution);

	assertSolution(solution, expected, numRows);

	delete[] solution;
}

static void test3x3Diagonal() {
	// [4  0  0] [x]   [8]
	// [0  3  0] [x] = [6]
	// [0  0  2] [x]   [4]
	int rowIndices[] = {0, 1, 2};
	int colIndices[] = {0, 1, 2};
	double values[] = {4.0, 3.0, 2.0};
	double rhs[] = {8.0, 6.0, 4.0};
	double expected[] = {2.0, 2.0, 2.0};

	runTest(3, 3, rowIndices, colIndices, values, rhs, expected);
}

static void test5x5Tridiagonal() {
	// [ 2 -1  0  0  0] [x]   [1]
	// [-1  2 -1  0  0] [x]   [0]
	// [ 0 -1  2 -1  0] [x] = [0]
	// [ 0  0 -1  2 -1] [x]   [0]
	// [ 0  0  0 -1  2] [x]   [1]
	int rowIndices[] = {0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4};
	int colIndices[] = {0, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4};
	double values[] = {2, -1, -1, 2, -1, -1, 2, -1, -1, 2, -1, -1, 2};
	double rhs[] = {1.0, 0.0, 0.0, 0.0, 1.0};
	double expected[] = {1.0, 1.0, 1.0, 1.0, 1.0};

	runTest(5, 13, rowIndices, colIndices, values, rhs, expected);
}

static void test5x5Sparse() {
	// [ 4  0 -1  0  0] [x]   [ 5.5]
	// [ 0  5  0 -2  0] [x]   [ 6.5]
	// [-1  0  6  0 -1] [x] = [ 0.5]
	// [ 0 -2  0  7  0] [x]   [16.0]
	// [ 0  0 -1  0  3] [x]   [ 2.5]
	int rowIndices[] = {0, 0, 1, 1, 2, 2, 2, 3, 3, 4, 4};
	int colIndices[] = {0, 2, 1, 3, 0, 2, 4, 1, 3, 2, 4};
	double values[] = {4, -1, 5, -2, -1, 6, -1, -2, 7, -1, 3};
	double rhs[] = {5.5, 6.5, 0.5, 16.0, 2.5};
	double expected[] = {1.5, 2.5, 0.5, 3.0, 1.0};

	runTest(5, 11, rowIndices, colIndices, values, rhs, expected);
}

// I generated a random 5x5 sparse matrix using Python scipy.sparse.random
static void test5x5Random() {
	// [ 0.000  0.375  0.000  0.916  0.000 ] [x]   [ 1.291 ]
	// [ 0.720  0.000  0.000  0.000  0.015 ] [x]   [ 0.735 ]
	// [ 0.000  0.000  0.468  0.823  0.000 ] [x] = [ 1.291 ]
	// [ 0.542  0.289  0.000  0.000  0.000 ] [x]   [ 0.831 ]
	// [ 0.000  0.112  0.634  0.000  0.000 ] [x]   [ 0.746 ]
	int numRows = 5;
	int numNonZeros = 10;

	int rowIndices[] = {0, 0, 1, 1, 2, 2, 3, 3, 4, 4};
	int colIndices[] = {1, 3, 0, 4, 2, 3, 0, 1, 1, 2};
	double values[] = {0.375, 0.916, 0.720, 0.015, 0.468,
										 0.823, 0.542, 0.289, 0.112, 0.634};
	double rhs[] = {1.291, 0.735, 1.291, 0.831, 0.746};
	double expected[] = {1.0, 1.0, 1.0, 1.0, 1.0};

	runTest(numRows, numNonZeros, rowIndices, colIndices, values, rhs, expected);
}

int main() {
	test3x3Diagonal();
	cout << "Test 3x3Diagonal passed.\n";

	test5x5Tridiagonal();
	cout << "Test 5x5Tridiagonal passed.\n";

	test5x5Sparse();
	cout << "Test 5x5Sparse passed.\n";

	test5x5Random();
	cout << "Test 5x5Random passed.\n";

	cout << "\nAll tests passed.\n";
	return 0;
}
