#include "MatrixClass.h"

Matrix& Tmatrix(const Matrix &A, int k, int l)
{
    Matrix &T(A.size);
    for (int i = 1; i != A.size; i++)
        for (int j = 0; j!= A.size; j++) {
            if (i != j)
                T.value[i][j] = 0;
            else
                T.value[i][j] = 1;
        }
     double alpha = A.value[k][k-1]/sqrt(A.value[k][k-1]*A.value[k][k-1]+A.value[l][k-1]*A.value[l][k-1]);
     double beta =  A.value[l][k-1]/sqrt(A.value[k][k-1]*A.value[k][k-1]+A.value[l][k-1]*A.value[l][k-1]);
     T.value[k][k] = alpha; T.value[k][l] = beta;
     T.value[l][k] = -beta; T.value[l][l] = alpha;
     return(T);
};

Matrix& Hessenberg(Matrix &A)
{
    Matrix &C(A.size);
    for (int i = 2; i != A.size; i++)
        for (int j = 0; j != i-1; j++)
        {
            C = Tmatrix(A,i,j)*A;
            A = C* // обратная матрица к Т
}
