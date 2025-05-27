#include <iostream>
#include <vector>
#include "matrix.h"

/* TO DO:
    fix up const operator overloads so we can call Matrices safely without altering them
    add consts to function parameters where possible
    make a CMAKE 
    add spacing in matrix indexing funcs for readability
    Decide between variant A or B for velocity projection
        - if variant A: ensure functions returns a struct
        - if variant B: ensure functions alters struct in place
        - in either case change all functions!!
    consider adding pressure matrix to your struct !!!!!!!!!
    consider a parameters.txt file for physical params
*/

/* HELPER FUNCS:
    assert that u and v have same size in funcs where we use them both
    assert that p and u/v have same size in funcs where we use them all
*/



// -----------------------------------------------------------------------------------------
//     PARAMETERS STRUCT
// -----------------------------------------------------------------------------------------
struct Velocity_2D 
{
    Matrix u;
    Matrix v;
};

struct FlowField2D
{
    Matrix u; // Velocity in the u direction
    Matrix v; // Velocity in the v direction
    Matrix p; // Pressure matrix
};

// -----------------------------------------------------------------------------------------
//     VELOCITY PROJECTION &
// -----------------------------------------------------------------------------------------

// // Variant A: return the whole struct by value
// Velocity_2D project_velocity(double rho, Matrix& u_in, Matrix& v_in, double dt, int dx, int dy, Matrix& p)
// {
//     Matrix u = u_in;
//     Matrix v = v_in;

//     int imax = p.getRows();
//     int jmax = p.getCols();

//     for (int i = 1; i < imax - 1; i++) {
//         for (int j = 1; j < jmax - 1; j++) {
//             u(i,j) -= dt * (1.0 / rho) * ( (p(i+1, j) - p(i-1, j)) / (2 * dx) );
//             v(i,j) -= dt * (1.0 / rho) * ( (p(i, j+1) - p(i, j-1)) / (2 * dy) );
//         }
//     }
//     return {u, v};
// }


// // Variant B: update an existing struct in place
// void project_velocity(double rho, double dt, int dx, int dy, Matrix& p, Velocity_2D& vel)
// {
//     Matrix& u = vel.u;
//     Matrix& v = vel.v;

//     const int imax = p.getRows();
//     const int jmax = p.getCols();

//     for (int i = 1; i < imax - 1; ++i) {
//         for (int j = 1; j < jmax - 1; ++j) {
//             u(i, j) -= dt * (1.0 / rho) * (p(i + 1, j) - p(i - 1, j)) / (2 * dx);
//             v(i, j) -= dt * (1.0 / rho) * (p(i, j + 1) - p(i, j - 1)) / (2 * dy);
//         }
//     }
// }


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
//     INTERMEDIATE VELOCITY 
// -----------------------------------------------------------------------------------------

// Variant A: return whole velocity struct by value
Velocity_2D calculate_intermediate_velocity(double nu, Matrix u, Matrix v, Matrix& u_old, Matrix& v_old, double dt, int dx, int dy)
{
    int imax = u.getRows();
    int jmax = u.getCols();

    for (int i = 1; i < imax - 1; ++i) {
        for (int j = 1; j < jmax - 1; ++j) {
            // Viscous (diffusive) term first
            u(i, j) = u_old(i, j) + dt * nu * ( (u_old(i + 1, j) + u_old(i - 1, j) - 2 * u_old(i, j)) / (dx * dx) + u_old(i, j + 1) + u_old(i, j - 1) - 2 * u_old(i, j) / (dy * dy));
            v(i, j) = v_old(i, j) + dt * nu * ( (v_old(i + 1, j) + v_old(i - 1, j) - 2 * v_old(i, j)) / (dx * dx) + v_old(i, j + 1) + v_old(i, j - 1) - 2 * v_old(i, j) / (dy * dy));

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
    return {u, v};
}


// Variant B_1: update existing velocity struct in place (uses u_old, v_old)
void calculate_intermediate_velocity(double nu, double dt, int dx, int dy, Velocity_2D& vel) 
{
    Matrix& u = vel.u;
    Matrix& v = vel.v;

    const int imax = u.getRows();
    const int jmax = u.getCols();


    for (int i = 1; i < imax - 1; ++i) {
        for (int j = 1; j < jmax - 1; ++j) {
            // Viscous (diffusive) term first
            u(i, j) = u_old(i, j) + dt * nu * ( (u_old(i + 1, j) + u_old(i - 1, j) - 2 * u_old(i, j)) / (dx * dx) + u_old(i, j + 1) + u_old(i, j - 1) - 2 * u_old(i, j) / (dy * dy));
            v(i, j) = v_old(i, j) + dt * nu * ( (v_old(i + 1, j) + v_old(i - 1, j) - 2 * v_old(i, j)) / (dx * dx) + v_old(i, j + 1) + v_old(i, j - 1) - 2 * v_old(i, j) / (dy * dy));

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


// Variant B_2: update existing velocity struct in place (uses u, v) 
// WRONG !!!!
void calculate_intermediate_velocity(double nu, double dt, int dx, int dy, Velocity_2D& vel) 
{
    Matrix& u = vel.u;
    Matrix& v = vel.v;

    const int imax = u.getRows();
    const int jmax = u.getCols();


    for (int i = 1; i < imax - 1; ++i) {
        for (int j = 1; j < jmax - 1; ++j) {
            // Viscous (diffusive) term first
            u(i, j) = u(i, j) + dt * nu * ( (u(i + 1, j) + u(i - 1, j) - 2 * u(i, j)) / (dx * dx) + u(i, j + 1) + u(i, j - 1) - 2 * u(i, j) / (dy * dy));
            v(i, j) = v(i, j) + dt * nu * ( (v(i + 1, j) + v(i - 1, j) - 2 * v(i, j)) / (dx * dx) + v(i, j + 1) + v(i, j - 1) - 2 * v(i, j) / (dy * dy));

            // Add the momentum (advection) terms using upwinding
            if (u(i, j) > 0.0) {
                u(i, j) -= dt * (u(i, j) * (u(i, j) - u(i - 1, j)) / dx);
                v(i, j) -= dt * (u(i, j) * (v(i, j) - v(i - 1, j)) / dx);
            }
            else {
                u(i, j) -= dt * (u(i, j) * (u(i + 1, j) - u(i, j)) / dx);
                v(i, j) -= dt * (u(i, j) * (v(i + 1, j) - v(i, j)) / dx);
            }

            if (v(i, j) > 0.0) {
                u(i, j) -= dt * (v(i, j) * (u(i, j) - u(i, j - 1)) / dy);
                v(i, j) -= dt * (v(i, j) * (v(i, j) - v(i, j - 1)) / dy);
            }
            else {
                u(i, j) -= dt * (v(i, j) * (u(i, j + 1) - u(i, j)) / dy);
                v(i, j) -= dt * (v(i, j) * (v(i, j + 1) - v(i, j)) / dy);
            }
        }
    }
}


// Variant B_3: update existing velocity struct in place and std::swap - WINNER WINNER CHICKEN DINNER
void calculate_intermediate_velocity(double nu, double dt, int dx, int dy, Velocity_2D& vel, Matrix& uScratch, Matrix& vScratch) 
{
    Matrix& u_old = vel.u;
    Matrix& v_old = vel.v;

    Matrix& u = uScratch;
    Matrix& v = vScratch;

    const int imax = u_old.getRows();
    const int jmax = u_old.getCols();


    for (int i = 1; i < imax - 1; ++i) {
        for (int j = 1; j < jmax - 1; ++j) {
            // Viscous (diffusive) term first
            u(i, j) = u_old(i, j) + dt * nu * ( (u_old(i + 1, j) + u_old(i - 1, j) - 2 * u_old(i, j)) / (dx * dx) + u_old(i, j + 1) + u_old(i, j - 1) - 2 * u_old(i, j) / (dy * dy));
            v(i, j) = v_old(i, j) + dt * nu * ( (v_old(i + 1, j) + v_old(i - 1, j) - 2 * v_old(i, j)) / (dx * dx) + v_old(i, j + 1) + v_old(i, j - 1) - 2 * v_old(i, j) / (dy * dy));

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
    std::swap(u, u_old);
    std::swap(v, v_old);
}



// Variant C: update existing velocity struct in place and std::swap with our updated struct - WINNER WINNER CHICKEN DINNER
void calculate_intermediate_velocity(double nu, double dt, int dx, int dy, FlowField2D& ff, Matrix& u_initial, Matrix& v_initial)
{
    std::swap(ff.u, u_initial);
    std::swap(ff.v, v_initial);

    Matrix& u = ff.u;
    Matrix& v = ff.v;

    const Matrix& u_old = u_initial;
    const Matrix& v_old = v_initial;

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
//     PRESSURE POISSON RIGHT HAND SIDE &
// -----------------------------------------------------------------------------------------

// Matrix calculate_ppm_RHS(double rho, Matrix u, Matrix v, double dt, int dx, int dy, Matrix p) 
// {
//     int imax = p.getRows();
//     int jmax = p.getCols();
//     Matrix RHS(imax,jmax, 0.0);  // ASSERT TO ENSURE THAT (imax, jmax) == X.dimensions!!!!!

//     for (int i = 1; i < imax - 1; i++) {
//         for (int j = 1; j < jmax - 1; j++) {
//             RHS(i,j) = rho * ((1.0 / dt) * ( (u(1+1, j) - u(i-1, j)) / (2 * dx) + (v(i, j+1) - v(i, j-1)) / (2 * dy))); 
//         }
//     }
//     return RHS;
// }

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
//     POISSON EQUATION SOLVER &
// -----------------------------------------------------------------------------------------

// Matrix pressure_poisson_jacobi(Matrix p, int dx, int dy, Matrix RHS, double rtol=1.0e-5, bool logs=false) 
// {
//     // code only works if dx == dy
//     double tol = 10.0 * rtol; 
//     int it = 0;
//     int imax = p.getRows();
//     int jmax = p.getCols();
//     Matrix p_old = p;
//     Matrix temp;

//     while (tol > rtol) {
//         it++;
//         temp = p_old;
//         p_old = p;
//         p = temp;

//         for (int i = 1; i < imax - 1; i++) {
//             for (int j = 1; j < jmax - 1; j++) {
//                 p(i,j) = 0.25 * (p_old(i+1, j) + p_old(i-1, j) + p_old(i, j+1) + p_old(i, j-1) - (dx * dx) * RHS(i, j));
//             }
//         }

//     // Apply zero pressure Neumann boundary conditions (dp/dy = 0)
//         // Apply BC's on left and right boundaries
//         for (int i = 0; i < imax; i++) {
//             p(i, 0) = p(i, 1);
//             p(i, jmax-1) = p(i, jmax-2);
//         }

//         // Apply BC's on top and bottom boundaries
//         for (int j = 0; j < jmax; j++) {
//             p(0, j) = p(1, j);
//             p(imax-1, j) = p(imax-2, j);
//         }

//         // Fix pressure level - choose arbitrary node to set p = 0
//         p(1,1) = 0.0;

//         // Relative change in pressure
//         tol = (p - p_old).Norm() / std::max(1.0e-10, p.Norm());    // Change Norm function to be static in Matrix.cpp
//     }

//     if (logs) {
//         std::cout << "Pressure solver iterations: " << it << std::endl;
//     }

//     return p;
// }


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
//     NAVIER STOKES SOLVER
// -----------------------------------------------------------------------------------------

void solve_NavierStokes(Velocity_2D& vel, Matrix& uScratch, Matrix& vScratch, Matrix& p, double rho, double nu, double dt, double t_end, int dx, int dy, double rtol=1.e-5, bool logs=false) 
{
    double t = 0.0;

    Matrix& u = vel.u;
    Matrix& v = vel.v;

    // Copy operation is fine here as it only happens once
    Matrix u_old(u);
    Matrix v_old(v);

    while (t < t_end) {
        t += dt;
        if (logs) std::cout << "Time = " << t << std::endl;

        // Calculate intermediate velocities
        calculate_intermediate_velocity(nu, dt, dx, dy, vel, uScratch, vScratch);

        // Calculate RHS for pressure poisson equation
        Matrix p_RHS = calculate_ppm_RHS(rho, u, v, dt, dx, dy, p);

        // Compute pressure
        p = pressure_poisson_jacobi()

        // double rho, Matrix u, Matrix v, double dt, int dx, int dy, Matrix p


    }
}

