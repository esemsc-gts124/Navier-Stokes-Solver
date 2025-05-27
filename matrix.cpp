#include "matrix.h"

using namespace std;

// YOU NEED TO FIX INDEX LOOPING MISMATCH FOR ROWS AND COLUMNS!!!!

// Constructor for any matrix if size and values are known
Matrix::Matrix(unsigned rowSize, unsigned colSize, double initial) {
    m_rowSize = rowSize;
    m_colSize = colSize;
    m_matrix.resize(rowSize);
    for (int i = 0; i < m_matrix.size(); i++)
        m_matrix[i].resize(colSize, initial);
}

// Default empty constructor for a matrix
Matrix::Matrix() {
    m_rowSize = 0;
    m_colSize = 0;
}


// Copy constructor
Matrix::Matrix(const Matrix &B) {
    this->m_colSize = B.getCols();
    this->m_rowSize = B.getRows();
    this->m_matrix = B.m_matrix;
}

// Destructor
Matrix::~Matrix() {}


// Matrix addition 
Matrix Matrix::operator+(Matrix &B) {
    Matrix sum(m_colSize, m_rowSize, 0.0); // ------- SWITCH AROUND ROW AND COL ARGS HERE!!
    unsigned i, j;
    for (i = 0; i < m_rowSize; i++) {
        for (j = 0; j < m_rowSize; j++) {
            sum(i,j) = this->m_matrix[i][j] + B(i,j);
        }
    }
    return sum;
}


// Matrix subtraction
Matrix Matrix::operator-(Matrix &B) {
    Matrix diff(m_colSize, m_rowSize, 0.0); // ------- SWITCH AROUND ROW AND COL ARGS HERE!!
    unsigned i, j;
    for (i = 0; i < m_rowSize; i++) {
        for (j = 0; j < m_colSize; j++) {
            diff(i,j) = this->m_matrix[i][j] - B(i,j);
        }
    }
    return diff;
}


// Matrix multiplication
Matrix Matrix::operator*(Matrix &B) {
    Matrix product(m_rowSize, B.getCols(), 0.0);
    if (m_colSize == B.getRows()) {
        unsigned i, j, k;
        double temp = 0.0;
        for (i = 0; i < m_rowSize; i++) {
            for (j = 0; j < m_colSize; j++) {
                temp = 0.0;
                for (k = 0; k < m_colSize; k++) {
                    temp += m_matrix[i][k] * B(k, j);
                }
                product(i,j) = temp;
            }
        }
        return product;
    }
    else {
        cout << "Error: row-column size mismatch between matrix A and B";
        Matrix error(1, 1, -99999);
        return error;
    }
}


// Scalar addition
Matrix Matrix::operator+(double scalar) {
    Matrix result(m_rowSize, m_colSize, 0.0);
    unsigned i, j;
    for (i = 0; i < m_rowSize; i++) {
        for (j = 0; j < m_colSize; j++) {
            result(i,j) = this->m_matrix[i][j] + scalar;
        }
    }
    return result;
}


// Scalar subtraction
Matrix Matrix::operator-(double scalar) {
    Matrix result(m_rowSize, m_colSize, 0.0);
    unsigned i, j;
    for (i = 0; i < m_rowSize; i++) {
        for (j = 0; j < m_colSize; j++) {
            result(i,j) = this->m_matrix[i][j] - scalar;
        }
    }
    return result;
}


// Scalar multiplication
Matrix Matrix::operator*(double scalar) {
    Matrix result(m_rowSize, m_colSize, 0.0);
    unsigned i, j;
    for (i = 0; i < m_rowSize; i++) {
        for (j = 0; j < m_colSize; j++) {
            result(i,j) = this->m_matrix[i][j] * scalar;
        }
    }
    return result;
}


// Scalar division
Matrix Matrix::operator/(double scalar) {
    Matrix result(m_rowSize, m_colSize, 0.0);
    unsigned i, j;
    for (i = 0; i < m_rowSize; i++) {
        for (j = 0; j < m_colSize; j++) {
            result(i,j) = this->m_matrix[i][j] / scalar; 
        }
    }
    return result;
}


// Return value of a given location when asked in form of A(x,y)
double& Matrix::operator()(const unsigned &rowNo, const unsigned &colNo) {
    return this->m_matrix[rowNo][colNo];
}


// Return row size
unsigned Matrix::getRows() const {
    return this->m_rowSize;
}


// Return col size
unsigned Matrix::getCols() const {
    return this->m_colSize;
}


// Transpose matrices
Matrix Matrix::transpose() {
    Matrix Transpose(m_colSize, m_rowSize, 0.0);
    unsigned i, j;
    for (i = 0; i < m_colSize; i++) {
        for (j = 0; j < m_rowSize; j++) {
            Transpose(i,j) = this->m_matrix[j][i];
        }
    }
    return Transpose;
}


// Frobenius norm of a matrix
double Matrix::Norm() {
    double sum = 0;
    double norm;
    unsigned i, j;
    for (i = 0; i < m_colSize; i++) {
        for (j = 0; j < m_rowSize; j++) {
            sum += this->m_matrix[i][j] * m_matrix[i][j];
        }
    }
    norm = sqrt(sum);
    return norm;
}


// Print the matrices 
void Matrix::print() const {
    cout << "Matrix: " << endl;
    unsigned i, j;
    for (i = 0; i < m_rowSize; i++) {
        for (j = 0; j < m_colSize; j++) {
            cout << "[" << m_matrix[i][j] << "]";
        }
        cout << endl;
    }
}