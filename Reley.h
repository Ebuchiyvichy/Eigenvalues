#pragma once
#include <fstream>
#include <vector>
#include <cmath>
#include "newton.h"

std::vector<double>	Newton_reley(int func, std::vector<double> x)
{
	std::vector<double>	xk(NBR);
	Matrix	J_(NBR);
	Matrix	R(NBR);
	Matrix	T(NBR);
	std::vector<double>	a(NBR);
	std::vector<double>	b(NBR);
	Matrix	J = Jacoby_matr(func, x);
	double	nor;

	R = J;
	T.onebyone();
	T.QR_find_x(R);
	J_.inverse_matrix(R, T);
	boaders(func, &a[0], &b[0]);
	do
	{
		iter++;
		if (iter > 30)
			return (x);
		xk = x;
		if (xk[0] < -a[0] || xk[0] > b[0] || xk[1] < -a[0] || xk[1] > b[0])
		{
	//		std::cout << "You are out of diap!" << std::endl;
			iter = 31;
			return (xk);
		}
		J = Jacoby_matr(func, xk);
		R = J;
		T.onebyone();
		T.QR_find_x(R);
		J_.inverse_matrix(R, T);
		for (int i = 0; i < NBR; i++)
			for (int j = 0; j < NBR; j++)
				if (abs(J_.value[i][j]) > 1.e+10)
				{
					iter = 31;
					return xk;
				}
		x = diff(xk, multi_vect(f(xk, func), J_));
	} while (norm(xk, x) > EPS);
//	for (int i = 0; i < NBR; i++)
//		std::cout << x[i] << '\t';
//	std::cout << std::endl;
	return x;
}

std::vector<double>	uniform_mesh(int func, int h1, char c)
{
	std::vector<double>	x(h1 + 1);
	double				a, b, h;

	boaders(func, &a, &b);
	if (c == 'a')
		h = (a + a) / h1;
	else
		h = (b + b) / h1;
	for (int i = 0; i <= h1; i++)
		x[i] = -a + h * i;
	return x;
}