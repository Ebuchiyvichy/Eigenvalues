#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

double	EPS = 10e-5;

class Matrix
{
public:
    std::vector<std::vector<double>> value;
    size_t size;

public:
    Matrix(size_t count)
    {
        size = count;
        for (int i = 0; i != size; i++)
        {
            value.push_back(std::vector<double>(count));
        }
    }
    ~Matrix(){

    }

    void init()
    {
        std::ifstream file;
        file.open("matrixexample.txt");

        for (int i = 0; i != size; i++)
            for (int j = 0; j != size; j++)
                file >> value[i][j];
        file.close();
    }
    void trunc()
    {
        double temp;
        for (int i = 0; i != size; i++)
        {
            for (int j = i+1; j != size; j++)
            {
                temp = value[i][j];
                value[i][j] = value[j][i];
                value[j][i] = temp;
            }
        }
    }
    void print()
    {
        for (int n = 0; n != size; n++)
        {
            std::cout << "\n";
            for (int j = 0; j != size; j++)
            {
                std::cout << std::setw(8) << value[n][j] << '\t';
            }
        }
        std::cout << "\n";
    }
    void onebyone()
    {
        for (int i = 0; i != size; i++)
            for (int j = 0; j != size; j++)
            {
                if (i == j)
                    value[i][i] = 1;
                else
                    value[i][j] = 0;
            }
    }
    void eigenvectors(const Matrix& A, const std::vector<double> lamb);
    friend void cpy(Matrix &A, Matrix& B)
    {
        for (int i = 0; i != A.size; i++)
            for (int j = 0; j != A.size; j++)
                B.value[i][j] = A.value[i][j];
    }
    void Relay (const Matrix& A, std::vector<double>& lamb);
    friend Matrix operator + (const Matrix &A, const Matrix &B)
    {
        Matrix C(A.size);
        if (A.size != B.size)
            std::cout << "Not correct size A and B";
        else
        {
            for (int i = 0; i != A.size; i++) {
                for (int j = 0; j != A.size; j++)
                    C.value[i][j] = A.value[i][j] + B.value[i][j];
            }
        }
        return C;
    }
    friend Matrix operator - (const Matrix &A, const Matrix &B)
    {
        Matrix C(A.size);
        if (A.size != B.size && C.size != B.size)
            std::cout << "Not correct size A and B";
        else
        {
            for (int i = 0; i != A.size; i++) {
                for (int j = 0; j != A.size; j++)
                    C.value[i][j] = A.value[i][j] - B.value[i][j];
            }
        }
        return C;
    }
    friend Matrix operator *(const Matrix &A, const Matrix &B)
    {
        Matrix C(A.size);
        for (int i = 0; i != A.size; i++)
            for (int j = 0; j != B.size; j++)
            {
                double sum = 0.0;
                for (int k = 0; k != A.size; k++)
                    sum += A.value[i][k]*B.value[k][j];
                C.value[i][j] = sum;
            }
        return C;
    }
    friend Matrix operator * (Matrix &A, double c)
    {
        for (int i = 0; i != A.size; i++) {
            for (int j = 0; j != A.size; j++)
                A.value[i][j] *=  c;
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
    friend std::vector<double> operator * (std::vector<double> &b, double c)
    {
        std::vector<double> a(b.size());
        for (int i = 0; i != b.size(); i++)
        {
            a[i] = b[i]*c;
        }
        return a;
    }
    friend double operator * (const std::vector<double>&a, const std::vector<double> &b)
    {
        double c = 0;
        for(int i = 0; i = a.size(); i++)
            c += a[i]*b[i];
        return c;
    }
    friend void delim(Matrix  &A, int k, std::vector<double>& x);
    friend void vych(Matrix &A, int k, std::vector<double>& x);
    friend void GaussRight(Matrix &A, std::vector<double>& x);
    friend void GaussLeft(const Matrix &A, std::vector <double>& x);
};
double skal (std::vector<double> a, std::vector <double> b)
{
    double sum = 0;
    for(int i = 0; i != a.size(); i++)
        sum += a[i]*b[i];
    return(sum);
}
//double pr(std::vector<double> a, std::vector<double> b)
//{
//    return (b * (skal(a,b)/skal(b,b)))
//}