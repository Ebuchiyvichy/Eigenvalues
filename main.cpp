#include "QR.h"

int main() {
    std::cout <<"Let's find eigenvalues" << std::endl;
    Matrix A(4);
    std::cout<< "Matrix A: "<< std::endl;
    A.init();

    A.print();
    A.value[1][0] = 5; A.value[0][1] = 5;
    std::vector<double> lamb(3);
    lamb = eigenvalues(A);
    std::cout<< "Eigenvalues: \n";
    for (int i = 0; i != A.size; i++)
       std::cout<< lamb[i] << '\t';
   return 0;
}
