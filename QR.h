//#include "Gauss.cpp"
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

bool buttomabs(const Matrix &A, int n)
{
	if (fabs(A.value[n][n - 1]) > pow(10, -4))
		return true;
	return false;
}

void qr(Matrix& A, Matrix& R, Matrix& T, int n)
{
	double c, s;
	double vspom;
	T.onebyone();
	R = A;

	for (int j = 0; j < n; j++)
	{
		for (int i = j + 1; i < n; i++)
		{
			if (fabs(R.value[i][j]) > EPS)
			{
				c = R.value[j][j] / sqrt((R.value[j][j]) * (R.value[j][j]) + (R.value[i][j]) * (R.value[i][j]));//вычисляем c12
				s = R.value[i][j] / sqrt((R.value[j][j]) * (R.value[j][j]) + (R.value[i][j]) * (R.value[i][j]));//вычисляем s12

				for (int k = 0; k < n; k++)
				{
					vspom = T.value[j][k];
					T.value[j][k] = T.value[j][k] * c + T.value[i][k] * s;
					T.value[i][k] = T.value[i][k] * c - vspom * s;
				}

				for (int k = 0; k < n; k++)//меняем матрицу А
				{
					vspom = R.value[j][k];
					R.value[j][k] = R.value[j][k] * c + R.value[i][k] * s;
					R.value[i][k] = R.value[i][k] * c - vspom * s;
				}

			}
		}
	}

	T.trunc();
}
//QR for all exc clear
/*
std::vector<double> eigenvalues(Matrix &A)
{
	Matrix T(A.size);
	Matrix R(A.size);
	Matrix Ak(A.size);
	Ak = A;
	std::vector<double> lambda(A.size);
	T.onebyone();
	int iterations = 0;
	Ak = Hessenberg(A);
	for (int n = A.size; n > 1; n--) {
		Ak.size = n;
		T.size = n;
		R.size = n;
		Matrix E(n);
		do {
			double delta = Ak.value[n - 1][n - 1];
			iterations++;
			E.onebyone();
			//Ak = Ak - E*delta;
			qr(Ak, R, T, n);
			Ak = R * T;
			E.onebyone();
			//Ak = Ak + E*delta;
		} while (buttomabs(Ak, n - 1));
		lambda[n - 1] = Ak.value[n - 1][n - 1];
	}
	std::cout << "Number of iterations: " << iterations << "\n";
	return lambda;
}*/

//void Matrix::eigenvectors(const Matrix& A, const std::vector<double> lamb)
//{
//    onebyone(); // matrix of our vectors
//    for (int i = 0; i != A.size; i++){
//        Matrix E(A.size); E.onebyone();
//        Matrix Ak(A.size); Ak = A - E*lamb[i];
//        std::vector<double> e_k(A.size);
//        do {
//            e_k = value[i];
//            GaussRight(Ak, value[i]);
//            GaussLeft(Ak, value[i]);
//            for (int j = 0; j != A.size; j++)
//                value[i][j] = value[i][j] / norm(value[i], A.size);
//        } while (fabs(fabs(skal(value[i], e_k)) -1) > EPS);
//    }
//
//
//}


bool buttomabs_clear(const Matrix &A, int n)
{
	for (int i = 1; i < n; i++)
		for (int j = 0; j < i; j++)
			if (fabs(A.value[i][j]) > pow(10, -4))
				return true;
	return false;
}
//QR for clear
std::vector<double> eigenvalues(Matrix &A)
{
	Matrix T(A.size);
	Matrix R(A.size);
	Matrix Ak(A.size);
	Ak = A;
	std::vector<double> lambda(A.size);
	T.onebyone();
	int iterations = 0;
	for (int n = A.size; n > 1; n--) {
		Ak.size = n;
		T.size = n;
		R.size = n;
		do {
			iterations++;
			qr(Ak, R, T, n);
			Ak = R * T;
		} while (buttomabs_clear(Ak, n - 1));
		lambda[n - 1] = Ak.value[n - 1][n - 1];
	}
	lambda[0] = Ak.value[0][0];
	std::cout << "Number of iterations: " << iterations << "\n";
	return lambda;
}