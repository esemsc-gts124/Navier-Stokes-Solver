#include <iostream>
#include <vector>
#include "matrix.h"

/* TO DO:
    add consts to function parameters where possible
    make a CMAKE 
    consider a parameters.txt file for physical params
*/

/* HELPER FUNCS:
    assert that u and v have same size in funcs where we use them both
    assert that p and u/v have same size in funcs where we use them all
*/



// -----------------------------------------------------------------------------------------
//     PARAMETERS STRUCT
// -----------------------------------------------------------------------------------------

struct FlowField2D
{
    Matrix u; // Velocity in the u direction
    Matrix v; // Velocity in the v direction
    Matrix p; // Pressure matrix
};



// -----------------------------------------------------------------------------------------
//     INTERMEDIATE VELOCITY 
// -----------------------------------------------------------------------------------------

// Variant C: update existing velocity struct in place and std::swap with our updated struct - WINNER WINNER CHICKEN DINNER
void calculate_intermediate_velocity(double nu, double dt, int dx, int dy, FlowField2D& ff, Matrix& u_initial, Matrix& v_initial)
{
    std::swap(ff.u, u_initial);
    std::swap(ff.v, v_initial);

    Matrix& u = ff.u;
    Matrix& v = ff.v;

    const Matrix& u_old = u_initial; // read-only
    const Matrix& v_old = v_initial; // read-only

    const int imax = u.getRows();
    const int jmax = u.getCols();

    for (int i = 1; i < imax - 1; ++i) {
        for (int j = 1; j < jmax - 1; ++j) {
            // Viscous (diffusive) term first
            u(i, j) = u_old(i, j) + dt * nu * ( ( u_old(i + 1, j) + u_old(i - 1, j) - 2 * u_old(i, j) ) / (dx * dx) + ( u_old(i, j + 1) + u_old(i, j - 1) - 2 * u_old(i, j) ) / (dy * dy) );
            v(i, j) = v_old(i, j) + dt * nu * ( ( v_old(i + 1, j) + v_old(i - 1, j) - 2 * v_old(i, j) ) / (dx * dx) + ( v_old(i, j + 1) + v_old(i, j - 1) - 2 * v_old(i, j) ) / (dy * dy) );

            // Add the momentum (advection) terms using upwinding
            if (u_old(i, j) > 0.0) {
                u(i, j) -= dt * (u_old(i, j) * (u_old(i, j) - u_old(i - 1, j)) / dx);
                v(i, j) -= dt * (u_old(i, j) * (v_old(i, j) - v_old(i - 1, j)) / dx);
            }
            else {
                u(i, j) -= dt * (u_old(i, j) * (u_old(i + 1, j) - u_old(i, j)) / dx);
                v(i, j) -= dt * (u_old(i, j) * (v_old(i + 1, j) - v_old(i, j)) / dx);
            }

            if (v_old(i, j) > 0.0) {
                u(i, j) -= dt * (v_old(i, j) * (u_old(i, j) - u_old(i, j - 1)) / dy);
                v(i, j) -= dt * (v_old(i, j) * (v_old(i, j) - v_old(i, j - 1)) / dy);
            }
            else {
                u(i, j) -= dt * (v_old(i, j) * (u_old(i, j + 1) - u_old(i, j)) / dy);
                v(i, j) -= dt * (v_old(i, j) * (v_old(i, j + 1) - v_old(i, j)) / dy);
            }
        }
    }
}



// -----------------------------------------------------------------------------------------
//     RIGHT HAND SIDE POISSON SOLVER
// -----------------------------------------------------------------------------------------

// Variant B: use our updated struct
Matrix calculate_ppm_RHS(double rho, double dt, int dx, int dy, FlowField2D ff)
{
    Matrix& u = ff.u;
    Matrix& v = ff.v;

    const int imax = u.getRows();
    const int jmax = u.getCols();
    Matrix RHS(imax, jmax, 0.0);

    for (int i = 1; i < imax - 1; ++i) {
        for (int j = 1; j < jmax - 1; ++j) {
            RHS(i, j) = rho * ( (1.0 / dt) * ( (u(i + 1, j) - u(i - 1, j)) / (2 * dx) + (v(i, j + 1) - v(i, j - 1)) / (2 * dy)));
        }
    }
    return RHS;
}



// -----------------------------------------------------------------------------------------
//     POISSON EQUATION SOLVER 
// -----------------------------------------------------------------------------------------

// Variant B: use our updated FF2D struct
void pressure_poisson_jacobi(int dx, int dy, Matrix& RHS, FlowField2D& ff, double rtol=1.e-5, bool logs=false)
{
    Matrix& p = ff.p;
    Matrix p_old = p;

    const int imax = p.getRows();
    const int jmax = p.getCols();

    int it = 0;
    double tol = 10.0 * rtol;

    while (tol > rtol) {
        ++it;

        for (int i = 1; i < imax - 1; ++i) {
            for (int j = 1; j < jmax - 1; ++j) {
                p(i, j) = 0.25 * (p_old(i + 1, j) + p_old(i - 1, j) + p_old(i, j + 1) + p_old(i, j - 1) - (dx * dx) * RHS(i, j));
            }
        }

        for (int i = 0; i < imax; ++i) {
            p(i, 0) = p(i, 1);
            p(i, jmax - 1) = p(i, jmax - 2);
        }

        for (int j = 0; j < jmax; ++j) {
            p(0, j) = p(1, j);
            p(imax - 1, j) = p(imax - 2, j);
        }

        p(1, 1) = 0.0;
        tol = (p - p_old).Norm() / std::max(1.e-10, p.Norm());
    
        std::swap(p, p_old);
    }

    if (logs) std::cout << "Pressure solver iterations: " << it << std::endl;
}



// -----------------------------------------------------------------------------------------
//     VELOCITY PROJECTION &
// -----------------------------------------------------------------------------------------

// Variant C: update an existing struct in place with our improved struct
void project_velocity(double rho, double dt, int dx, int dy, FlowField2D& ff)
{
    Matrix& u = ff.u;
    Matrix& v = ff.v;

    const int imax = u.getRows(); 
    const int jmax = u.getCols();

    for (int i = 1; i < imax - 1; ++i) {
        for (int j = 1; j < jmax - 1; ++j) {
            u(i, j) -= dt * (1.0 / rho) * ( p(i + 1, j) - p(i - 1, j)) / (2 * dy);
            v(i, j) -= dt * (1.0 / rho) * ( p(i, j + 1) - p(i, j - 1)) / (2 * dy);
        }
    }
}



// -----------------------------------------------------------------------------------------
//     NAVIER STOKES SOLVER
// -----------------------------------------------------------------------------------------

// void solve_NavierStokes(FlowField2D& ff, double rho, double nu, double dt, double t_end, int dx, int dy, double rtol=1.e-5, bool logs=false)
// {
//     Matrix u_initial(ff.u);
//     Matrix v_initial(ff.v);

//     Matrix p_RHS;
    
//     double t = 0.0;

//     while (t < t_end) {
//         t += dt;
//         if (logs) std::cout << "Time =  " << t << std::endl;

//         calculate_intermediate_velocity();
//         p_RHS = calculate_ppm_RHS();
//         pressure_poisson_jacobi();
//         project_velocity


//     }
// } 


