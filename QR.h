#pragma once
#include "MatrixClass.h"


Matrix Tmatrix(const Matrix &A, int k, int l)
{
	Matrix T(A.size);
	for (int i = 0; i != A.size; i++)
		for (int j = 0; j != A.size; j++) {
			if (i != j)
				T.value[i][j] = 0;
			else
				T.value[i][j] = 1;
		}
	double alpha = A.value[k][k - 1] / sqrt(A.value[k][k - 1] * A.value[k][k - 1] + A.value[l][k - 1] * A.value[l][k - 1]);
	double beta = A.value[l][k - 1] / sqrt(A.value[k][k - 1] * A.value[k][k - 1] + A.value[l][k - 1] * A.value[l][k - 1]);
	T.value[k][k] = alpha; T.value[k][l] = beta;
	T.value[l][k] = -beta; T.value[l][l] = alpha;
	return(T);
};

Matrix Hessenberg(Matrix &A)
{
	Matrix C(A.size);
	Matrix T(A.size);
	for (int j = 0; j != A.size - 2; j++)
		for (int i = j + 2; i != A.size; i++) {
			T = Tmatrix(A, j + 1, i);
			C = T * A;
			T.trunc();
			A = C * T;
		}
	return A;
}

void Matrix::QR(Matrix &A)
{
	double	c;
	double	s;
	double	tmp;

	for (int k = 0; k < size; k++)
	{
		for (int i = k + 1; i < size; i++)
		{
			if (fabs(A.value[k][i]) > 10e-8)
			{
				c = A.value[k][k] / sqrt(A.value[k][k] * A.value[k][k] + A.value[i][k] * A.value[i][k]);
				s = A.value[i][k] / sqrt(A.value[k][k] * A.value[k][k] + A.value[i][k] * A.value[i][k]);
				for (int j = 0; j <= i; j++) // change T-matrix
				{
					tmp = value[k][j];
					value[k][j] = value[k][j] * c + value[i][j] * s;
					value[i][j] = c * value[i][j] - s * tmp;
				}
				for (int j = k; j < size; j++) // change A-matrix
				{
					tmp = A.value[k][j];
					A.value[k][j] = c * A.value[k][j] + s * A.value[i][j];
					A.value[i][j] = c * A.value[i][j] - s * tmp;
				}
			}
		}
	}
	//    std::cout << "QR method matrix R:" << std::endl;
	//    A.print();
	//
	//    std::cout << "QR method matrix T:" << std::endl;
	//    print();
	//		trunc();
}

bool buttomabs(const Matrix &A)
{
	//    for(int i = 0; i != A.size; i++)
	//    {
	//        double sum = 0;
	//        for (int j = 0; j != i; j++)
	//            sum += fabs(A.value[i][j]);
	//        std::cout << sum <<"\n";
	//        if (sum >  pow(10,-4))
	//            return true;
	//    }
	//    return false;
	for (int i = 1; i != A.size; i++)
		if (fabs(A.value[i][i - 1]) > pow(10, -4))
			return true;
	return false;
}

std::vector<double> eigenvalues(Matrix &A)
{
	Matrix T(A.size);
	Matrix R(A.size);
	Matrix Ak(A.size);
	Ak = A;
	std::vector<double> lambda(A.size);
	T.onebyone();
	int iterations = 0;
	R = Ak;
	T.QR(R);
	T.trunc();
	Ak = R * T;
	//	Ak = Hessenberg(A);
	do {
		iterations++;
		R = Ak;
		T.QR(R);
		Ak = R * T;
	} while (buttomabs(Ak));
	for (int i = 0; i != A.size; i++)
		lambda[i] = Ak.value[i][i];
	Ak.print();
	//אכדמנטעל סמ
	std::cout << "Number of iterations: " << iterations << "\n";
	return lambda;
}