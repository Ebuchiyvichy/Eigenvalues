#include "QR.h"

int main() {
	std::cout << "Let's find eigen values" << std::endl;
	Matrix A(4);
	std::cout << "Matrix A: " << std::endl;
	A.init();
	A.print();


	std::vector<double> lamb(A.size);
	lamb = eigenvalues(A);
	std::cout << "Eigen values: \n";
	for (int i = 0; i != A.size; i++)
		std::cout << lamb[i] << '\t';

	system ("pause");
	return 0;
}