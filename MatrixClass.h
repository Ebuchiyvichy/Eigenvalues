#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

double	EPS = 10e-8;

class Matrix:
{
public:
    std::vector<std::vector<double>> value;
    size_t size;

public:
    Matrix(size_t count)
    {
        size = count;
        value.resize(count);
        for (int i = 0; i != size; i++)
        {
            value.push_back(std::vector<double>(count));
        }
        std::cout << "I can construct!\n";
    }
    ~Matrix(){

    }

    void init()
    {
        std::ifstream file;
        file.open("matrixexample.txt");

        for (int i = 0; i != size; i++)
            for (int j = 0; j <= size; j++)
            {
                double numd;
                file >>  numb;
                value[i][j] = numb;
            }
        file.close();
    }

    Matrix& operator + (const Matrix &A, const Matrix &B)
    {
        if (A.size != B.size && size != B.size)
            std::cout << "Not correct size A and B";
        else
        {
            for (int i = 0; i != A.size; i++) {
                for (int j = 0; j != A.size; j++)
                    value[i][j] = A.value[i][j] + B.value[i][j];
            }
        }
        return *this;
    }
    Matrix& operator - (const Matrix &A, const Matrix &B)
    {
        if (A.size != B.size && size != B.size)
            std::cout << "Not correct size A and B";
        else
        {
            for (int i = 0; i != A.size; i++) {
                for (int j = 0; j != A.size; j++)
                    value[i][j] = A.value[i][j] - B.value[i][j];
            }
        }
        return *this;
    }
    friend Matrix& operator *(const Matrix &A, const Matrix &B)
    {
        for (int i = 0; i != A.size; i++)
            for (int j = 0; j != B.size; j++)
            {
                double sum = 0;
                for (int k = 0; k != A.size; k++)
                    sum += A.value[j][k]*B.value[k][j];
                value[i][j] = sum;
            }
        return *this;
    }
    friend Matrix& operator * (Matrix &A, double c)
    {
        for (int i = 0; i != A.size; i++) {
            for (int j = 0; j != A.size; j++)
                A.value[i][j] = A.value[i][j] * c;
        }
        return A;
    }
    friend std::vector<double> operator * (const Matrix &A, const std::vector<double> &b) {
        std::vector<double> c(b.size());
        if (A.size != b.size())
            std::cout << "Wrong size in vector or matrix!\n";
        else {
            for (int i = 0; i != A.size; i++)
                for (int j = 0; j != A.size; j++) {
                    c[i] += A.value[i][j] * b[j];
                }
        }
        return c;
    }
};