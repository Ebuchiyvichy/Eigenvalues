#pragma once
#include <fstream>
#include <vector>
#include <cmath>
#include "MatrixClass.h"

double			EPS = 10.e-6;
static int		HBISEC = 7;
int				TEST = -2;
int	const		NBR = 2;	
int const		N = 300;	//количество разбиений для бассейна Ньютона
static int	iter;


typedef struct local
{
	int			nbr;
	double 		**intervals;
};

double	f(double x, int func)
{
	if (func == 1)
		return ((x - 0.1)*(x - 0.22)*(x - 0.55)*(x - 0.7)*(x - 0.75));
	else if (func == 2)
		return sqrt(x + 1) - 1;
	else if (func == 3)
		return 35 * pow(x, 3) - 67 * pow(x, 2) - 3 * x + 3;
	else
		return (0);
}

std::vector<double>	f(std::vector<double> x, int func)
{
	if (func == -1)
	{
		std::vector<double>	y(2);
		y[0] = x[0] * x[0] - x[1] * x[1] - 15;
		y[1] = x[0] * x[1] + 4;
		return y;
	}
	else if (func == -2)
	{
		std::vector<double>	y(2);
		y[0] = x[0] * x[0] + x[1] * x[1] + x[0] + x[1] - 8;
		y[1] = x[0] * x[0] + x[1] * x[1] + x[0] * x[1] - 7;
		return y;
	}
}

void	boaders(int func, double *a, double *b)
{
	if (func == 1)
	{
		*a = 0.0;
		*b = 1.0;
	}
	else if (func == 2)
	{
		*a = -1.0;
		*b = 10.0;
	}
	else if (func == 3)
	{
		*a = 0.0;
		*b = 1.0;
	}
	else if (func == -1)
	{
		*a = 10;
		*b = 10;
	}
	else if (func == -2)
	{
		*a = 10;
		*b = 10;
	}
}

void Matrix::QR_find_x(Matrix& A)
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
}

std::vector<double>	cpy_vector(std::vector<double> tmp, std::vector<double> x, int size)
{

	for (int i = 0; i < size; i++)
		tmp[i] = x[i];
	return (tmp);
}

std::vector<double>	find_x(const Matrix& A, std::vector<double> b_)
{
	std::vector<double>	x(A.size);

	for (int i = A.size - 1; i >= 0; i--)
	{
		for (int j = i + 1; j < A.size; j++)
			b_[i] = b_[i] - A.value[i][j] * b_[j];
		b_[i] = b_[i] / A.value[i][i];
	}
	x = cpy_vector(x, b_, A.size);
	return (x);
}

void Matrix::inverse_matrix(const Matrix& R, Matrix& T)
{
	//with Hausse
	std::vector<double> x(T.size);
	std::vector<double> y(T.size);
	Matrix E(R.size);

	E.onebyone();
	for (int i = 0; i < T.size; i++)
	{
		for (int j = 0; j < T.size; j++)
		{
			if (j == i)
				y[j] = 1;
			else
				y[j] = 0;
		}
		x = find_x(R, multi_vect(y, T));
		for (int j = 0; j < T.size; j++)
		{
			value[i][j] = x[j];
		}
	}
	trunc();
}

std::vector<double>	multi_vect(std::vector<double> I, const Matrix& T)
{
	std::vector<double> tmp(NBR);

	for (int j = 0; j < T.size; j++)
		for (int i = 0; i < T.size; i++)
			tmp[j] += T.value[j][i] * I[i];
	return (tmp);
}

std::vector<double>	diff(std::vector<double> a, std::vector<double> b)
{
	std::vector<double>	tmp(NBR);

	for (int i = 0; i < NBR; i++)
		tmp[i] = a[i] - b[i];
	return (tmp);
}

double	norm(std::vector<double> a, std::vector<double> b)
{
	double	max;
	max = fabs(a[0] - b[0]);
	for (int i = 1; i < NBR; i++)
		if (fabs(a[i] - b[i]) > max)
			max = fabs(a[i] - b[i]);
	return max;
}