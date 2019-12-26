#pragma once
#include <fstream>
#include <vector>
#include <cmath>
#include "functions.h"

local	localisation_func(int func, int nbrh)
{
	local	l;
	double	a, b;
	double	h;
	double	x;
	
	l.nbr = 0;
	boaders(func, &a, &b);
	h = fabs(b - a) / HBISEC;
	l.intervals = new double*[nbrh];	/*создаю лишние ячейки*/
	for (int i = 0; i <= nbrh; i++)
		l.intervals[i] = new double[2];
	for (int i = 0; i < nbrh; i++)
	{
		x = a + h * i;
		if (f(x, func) * f(x + h, func) <= 0)
		{
			l.intervals[l.nbr][0] = x;
			l.intervals[l.nbr][1] = x + h;
			l.nbr++;
		}
	}
	return l;
}

double	Bisection(int func, local l, int n)
{
	double	x, a, b;
	int		i = 1;

	a = l.intervals[n][0];
	b = l.intervals[n][1];
	x = (a + b) / 2;
	while (fabs(b - a) > EPS)
	{
		i++;
		if (f(a, func) * f(x, func) <= 0)
			b = x;
		else
			a = x;
		x = (a + b) / 2;
	}
	std::cout << "Nbr of iteration in bisection: " << i << std::endl;
	return x;
}

double	Hords(int func, local l, int n)
{
	double	x1, x2, a, b;
	int		i = 1;

	a = l.intervals[n][0];
	b = l.intervals[n][1];

	x1 = (f(a, func) * b - f(b, func) * a) / (f(a, func) - f(b, func));
	do
	{
		i++;
		x2 = x1;
		x1 = x2 - f(x2, func)*(b - x2) / (f(b, func) - f(x2, func));
	} while (fabs(x2 - x1) > EPS);
	std::cout << "Nbr of iteration in hords: " << i << std::endl;
	return (x1);
}

double	D(int func, double x)
{
	return (f(x + EPS, func) - f(x, func)) / EPS;
}

double	D_cen(int func, double x)
{
	return (f(x + EPS, func) - f(x - EPS, func)) / (2 * EPS);
}

double	Newton(int func, local l, int n)
{
	double	x1, x2, a, b;
	int		i = 1;

	a = l.intervals[n][0];
	b = l.intervals[n][1];

	x2 = (f(a, func) * b - f(b, func) * a) / (f(a, func) - f(b, func));
	do
	{
		i++;
		x1 = x2;
		if (D_cen(func, x1) < EPS)
		{
			std::cout << "Trouble" << std::endl;
			return (Hords(func, l, n));
		}
		x2 = x1 - f(x1, func) / D_cen(func, x1);
	} while (fabs(x1 - x2) > EPS);
	std::cout << "Nbr of iteration in Newton: " << i << std::endl;
	return (x2);
}
