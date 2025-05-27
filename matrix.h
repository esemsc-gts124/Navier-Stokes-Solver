#ifndef __EE_242_Project_2__matrix__
#define __EE_242_Project_2__matrix__

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <cmath>

using std::vector;
using std::tuple;

class Matrix {
    private:
        unsigned m_rowSize;
        unsigned m_colSize;
        vector<vector<double> > m_matrix;

    public:
        Matrix(unsigned, unsigned, double);
        // Matrix(const char *); // This one is for reading from a file
        Matrix(const Matrix &);
        Matrix();
        ~Matrix();

        // Matrix Operations
        Matrix operator+(Matrix &);
        Matrix operator-(Matrix &);
        Matrix operator*(Matrix &);
        Matrix transpose();
        double Norm();

        // Scalar Operations
        Matrix operator+(double);
        Matrix operator-(double);
        Matrix operator*(double);
        Matrix operator/(double);

        // Aestheic Methods
        double& operator()(const unsigned &, const unsigned &);
        void print() const;
        unsigned getRows() const;
        unsigned getCols() const;

        // Power Iteration
        tuple<Matrix, double, int> powerIter(unsigned, double);

        // Deflation
        Matrix deflation(Matrix &, double&);
};
#endif