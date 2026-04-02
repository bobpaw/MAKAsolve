#include <MAKAsolve/Solver.h>
#include <iostream>
#include <cmath>

using namespace std;

int main() {
  int numRows = 3;
  int numNonZeros = 3;

  // COO format
  int rowIndices[] = {0, 1, 2};
  int colIndices[] = {0, 1, 2};
  double values[] = {4.0, 3.0, 2.0};

  double rhs[] = {8.0, 6.0, 4.0};
  double solution[] = {0.0, 0.0, 0.0};

  maka::Solver s(numRows, numNonZeros, rowIndices, colIndices, values);
  s.solve(rhs, solution);

  cout << "Solution = {" << solution[0] << ", " << solution[1] << ", " << solution[2] << "}\n";
	cout << "Expected = {2, 2, 2}\n";

  return 0;
}