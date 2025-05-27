#include <iostream>
#include <vector>
#include "matrix.h"

/* TO DO:
    fix up const operator overloads so we can call Matrices safely without altering them
    sort out RHS matrix in ppm_RHS func
    make a CMAKE 
    add spacing in matrix indexing funcs for readability
    Decide between variant A or B for velocity projection
    consider adding pressure matrix to your struct
*/



struct Velocity_2D 
{
    Matrix u;
    Matrix v;
};

// Variant A: return the whole struct by value
Velocity_2D project_velocity(double rho, Matrix& u_in, Matrix& v_in, double dt, int dx, int dy, Matrix& p)
{
    Matrix u = u_in;
    Matrix v = v_in;

    int imax = p.getRows();
    int jmax = p.getCols();

    for (int i = 1; i < imax - 1; i++) {
        for (int j = 1; j < jmax - 1; j++) {
            u(i,j) -= dt * (1.0 / rho) * ( (p(i+1, j) - p(i-1, j)) / (2 * dx) );
            v(i,j) -= dt * (1.0 / rho) * ( (p(i, j+1) - p(i, j-1)) / (2 * dy) );
        }
    }
    return {u, v}
}


// Variant B: update an existing struct in place
void project_velocity(double rho, double dt, int dx, int dy, Matrix& p, Velocity_2D& vel)
{
    Matrix& u = vel.u;
    Matrix& v = vel.v;

    const int imax = p.getRows();
    const int jmax = p.getCols();

    for (int i = 1; i < imax - 1; ++i) {
        for (int j = 1; j < jmax - 1; ++j) {
            u(i, j) -= dt * (1.0 / rho) * (p(i + 1, j) - p(i - 1, j)) / (2 * dx);
            v(i, j) -= dt * (1.0 / rho) * (p(i, j + 1) - p(i, j - 1)) / (2 * dy);
        }
    }
}




Matrix calculate_ppm_RHS(double rho, Matrix u, Matrix v, double dt, int dx, int dy, Matrix p) 
{
    int imax = p.getRows();
    int jmax = p.getCols();
    Matrix RHS = (?,?, 0.0);  // I guess RHS size is determined by the size of your domain

    for (int i = 1; i < imax - 1; i++) {
        for (int j = 1; j < jmax - 1; j++) {
            RHS(i,j) = rho * ((1.0 / dt) * ( (u(1+1, j) - u(i-1, j)) / (2 * dx) + (v(i, j+1) - v(i, j-1)) / (2 * dy))); 
        }
    }
    return RHS;
}








Matrix pressure_poisson_jacobi(Matrix p, int dx, int dy, Matrix RHS, double rtol=1.0e-5, bool logs=false) 
{
    // code only works if dx == dy
    double tol = 10.0 * rtol; 
    int it = 0;
    int imax = p.getRows();
    int jmax = p.getCols();
    Matrix p_old = p;
    Matrix temp;

    while (tol > rtol) {
        it++;
        temp = p_old;
        p_old = p;
        p = temp;

        for (int i = 1; i < imax - 1; i++) {
            for (int j = 1; j < jmax - 1; j++) {
                p(i,j) = 0.25 * (p_old(i+1, j) + p_old(i-1, j) + p_old(i, j+1) + p_old(i, j-1) - (dx * dx) * RHS(i, j));
            }
        }

    // Apply zero pressure Neumann boundary conditions (dp/dy = 0)
        // Apply BC's on left and right boundaries
        for (int i = 0; i < imax; i++) {
            p(i, 0) = p(i, 1);
            p(i, jmax-1) = p(i, jmax-2);
        }

        // Apply BC's on top and bottom boundaries
        for (int j = 0; j < jmax; j++) {
            p(0, j) = p(1, j);
            p(imax-1, j) = p(imax-2, j);
        }

        // Fix pressure level - choose arbitrary node to set p = 0
        p(1,1) = 0.0;

        // Relative change in pressure
        tol = (p - p_old).Norm() / std::max(1.0e-10, p.Norm());    // Change Norm function to be static in Matrix.cpp
    }

    if (logs) {
        std::cout << "Pressure solver iterations: " << it << std::endl;
    }

    return p;
}




int main() {
    Matrix A(4,5,1.0);
    A.print();
    A(1,1) = 0.0;
    A.print();
}

