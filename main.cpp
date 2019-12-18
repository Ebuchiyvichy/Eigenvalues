#include "QR.h"

int main() {
    std::cout <<"Let's find eigenvalues" << std::endl;
    Matrix A(4);
    std::cout<< "Matrix A: "<< std::endl;
    A.init();
    A.print();
    std::vector<double> lamb(A.size);
    lamb = eigenvalues(A);
    std::cout<< "Eigenvalues: \n";
    for (int i = 0; i != A.size; i++)
       std::cout<< lamb[i] << '\t';
    std::cout << std::endl;
    A.init();
    std::cout << " Reverse Iterations:\n";
    Matrix eigenvect(4);
    eigenvect.eigenvectors(A, lamb);
    eigenvect * 0.5;
    for (int i = 0; i != A.size; i++)
        lamb[i] = 0;
    eigenvect.Relay(A,lamb);
    eigenvect.print();

   return 0;
}
