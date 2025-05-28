#include <iostream>
#include <vector>
#include <cmath>
#include "matrix.h"
#include <iomanip>
#include <algorithm>
#include <chrono>

/* TO DO:
    TEST OUT SWAPPING OPERATIONS INSTEAD OF COPY IN NS_SOLVER!! TIME IT!!
    add consts to function parameters where possible
    make a CMAKE 
    consider a parameters.txt file for physical params
    change dx and dy to doubles to avoid integer division
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
void calculate_intermediate_velocity(double nu, double dt, double dx, double dy, FlowField2D& ff, Matrix& u_initial, Matrix& v_initial)
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
            u(i, j) = u_old(i, j) + dt * nu * ( ( u_old(i + 1, j) + u_old(i - 1, j) - 2 * u_old(i, j) ) / (dx * dx) + 
            ( (u_old(i, j + 1) + u_old(i, j - 1) - 2 * u_old(i, j)) ) / (dy * dy) );

            v(i, j) = v_old(i, j) + dt * nu * ( ( v_old(i + 1, j) + v_old(i - 1, j) - 2 * v_old(i, j) ) / (dx * dx) + 
            ( (v_old(i, j + 1) + v_old(i, j - 1) - 2 * v_old(i, j)) ) / (dy * dy) );

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
Matrix calculate_ppm_RHS(double rho, double dt, double dx, double dy, const FlowField2D& ff)
{
    const Matrix& u = ff.u;
    const Matrix& v = ff.v;

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
void pressure_poisson_jacobi(double dx, double dy, Matrix& RHS, FlowField2D& ff, double rtol=1.e-5, bool logs=false)
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

    std::swap(p, p_old);

    if (logs) std::cout << "Pressure solver iterations: " << it << std::endl;
}



// -----------------------------------------------------------------------------------------
//     VELOCITY PROJECTION 
// -----------------------------------------------------------------------------------------

// Variant C: update an existing struct in place with our improved struct
void project_velocity(double rho, double dt, double dx, double dy, FlowField2D& ff)
{
    Matrix& u = ff.u;
    Matrix& v = ff.v;
    const Matrix& p = ff.p;

    const int imax = u.getRows(); 
    const int jmax = u.getCols();

    for (int i = 1; i < imax - 1; ++i) {
        for (int j = 1; j < jmax - 1; ++j) {
            u(i, j) -= dt * (1.0 / rho) * ( p(i + 1, j) - p(i - 1, j)) / (2 * dx);
            v(i, j) -= dt * (1.0 / rho) * ( p(i, j + 1) - p(i, j - 1)) / (2 * dy);
        }
    }
}



// -----------------------------------------------------------------------------------------
//     REAPPLY BOUNDARY CONDITIONS
// -----------------------------------------------------------------------------------------
inline void applyBoundaryConditions(FlowField2D& ff, unsigned Nx, unsigned Ny)
{
    // for (unsigned j = 0; j < Ny; ++j)
    // {
    //     /* left wall */
    //     ff.u(0,   j)  = 0.0;
    //     ff.v(0,   j)  = 0.0;

    //     /* right wall */
    //     ff.u(Nx-1,j)  = 0.0;
    //     ff.v(Nx-1,j)  = 0.0;
    // }

    // for (unsigned i = 0; i < Nx; ++i)
    // {
    //     /* lid: u = 1, v = 0 */
    //     ff.u(i, Ny-1) = 1.0;
    //     ff.v(i, Ny-1) = 0.0;

    //     /* bottom wall */
    //     ff.u(i, 0)    = 0.0;
    //     ff.v(i, 0)    = 0.0;
    // }
    // /* left & right walls ------------------------------------------------ */
    // for (unsigned j = 0; j < Ny; ++j) {
    //     ff.u(0,   j) = ff.v(0,   j) = 0.0;
    //     ff.u(Nx-1,j) = ff.v(Nx-1,j) = 0.0;
    // }

    // /* bottom wall ------------------------------------------------------- */
    // for (unsigned i = 0; i < Nx; ++i) {
    //     ff.u(i, 0) = ff.v(i, 0) = 0.0;
    // }

    // /* lid (must come last) --------------------------------------------- */
    // for (unsigned i = 0; i < Nx; ++i) {
    //     ff.u(i, Ny-1) = 1.0;
    //     ff.v(i, Ny-1) = 0.0;
    // }

    // Left and right walls
    for (int j = 0; j < Ny; ++j) {
        ff.u(0, j) = 0.0;
        ff.v(0, j) = 0.0;

        ff.u(Nx-1, j) = 0.0;
        ff.v(Nx-1, j) = 0.0;
    }

    // Bottom wall
    for (int i = 0; i < Nx; ++i) {
        ff.u(i, 0) = 0.0;
        ff.v(i, 0) = 0.0;
    }

    // Lid (must be the last to be set or else we risk overwriting the corner values)
    for (int i = 0; i < Nx; ++i) {
        ff.u(i, Ny-1) = 1.0;
        ff.v(i, Ny-1) = 0.0;
    }
}



// -----------------------------------------------------------------------------------------
//     NAVIER STOKES SOLVER
// -----------------------------------------------------------------------------------------

// void solve_NavierStokes(FlowField2D& ff, double rho, double nu, double dt, double t_end, double dx, double dy, double rtol=1.e-5, bool logs=false)
// {
//     Matrix u_initial(ff.u);
//     Matrix v_initial(ff.v);

//     Matrix p_RHS;
//     double tolu;
//     double tolv;
    
//     double t = 0.0;

//     while (t < t_end) {
//         t += dt;
//         if (logs) std::cout << "Time =  " << t << std::endl;

//         calculate_intermediate_velocity(nu, dt, dx, dy, ff, u_initial, v_initial);
//         p_RHS = calculate_ppm_RHS(rho, dt, dx, dy, ff);
//         pressure_poisson_jacobi(dx, dy, p_RHS, ff, rtol, logs);
//         project_velocity(rho, dt, dx, dy, ff);

//         if (logs) {
//             std::cout << "Norm(u): " << ff.u.Norm() << std::setprecision(8) << "    Norm(v): " << ff.v.Norm() << std::setprecision(8) << std::endl;
//             // std::cout << "Courant number: " << std::max(sqrt)
//         }

//         tolu = (ff.u - u_initial).Norm() / std::max(1.0e-10, ff.u.Norm());
//         tolv = (ff.v - v_initial).Norm() / std::max(1.0e-10, ff.v.Norm());

//         if (tolu < rtol && tolv < rtol) break;

//         std::swap(u_initial, ff.u);
//         std::swap(v_initial, ff.v);
//     }
// } 

void solve_NavierStokes(FlowField2D& ff, double rho, double nu, double dt, double t_end, double dx, double dy, double Nx, double Ny, double rtol=1.e-5, bool logs=false)
{
    Matrix u_initial(ff.u);
    Matrix v_initial(ff.v);

    Matrix p_RHS;
    double tolu;
    double tolv;
    
    double t = 0.0;

    while (t < t_end) {
        t += dt;
        if (logs) std::cout << "Time =  " << t << std::endl;

        calculate_intermediate_velocity(nu, dt, dx, dy, ff, u_initial, v_initial);
        p_RHS = calculate_ppm_RHS(rho, dt, dx, dy, ff);
        pressure_poisson_jacobi(dx, dy, p_RHS, ff, rtol, logs);
        project_velocity(rho, dt, dx, dy, ff);
        applyBoundaryConditions(ff, Nx, Ny);

        if (logs) {
            std::cout << "Norm(u): " << ff.u.Norm() << std::setprecision(8) << "    Norm(v): " << ff.v.Norm() << std::setprecision(8) << std::endl;
            // std::cout << "Courant number: " << std::max(sqrt)
        }

        tolu = (ff.u - u_initial).Norm() / std::max(1.0e-10, ff.u.Norm());
        tolv = (ff.v - v_initial).Norm() / std::max(1.0e-10, ff.v.Norm());

        if (tolu < rtol && tolv < rtol) break;

        // std::swap(u_initial, ff.u);
        // std::swap(v_initial, ff.v);

        u_initial = ff.u;
        v_initial = ff.v;
    }
} 



int main() 
{
    // Physical constants
    double rho = 1.0;
    double nu = 1.0 / 10.0;

    // Spatial mesh domain
    double Lx = 1;
    double Ly = 1;

    // Number of grid points including boundary nodes
    double Nx = 51;
    double Ny = 51;

    // Mesh spacing
    double dx = Lx / (Nx - 1);
    double dy = Ly / (Ny - 1);

    // INSERT X AND Y GRIDS HERE
    Matrix X(Nx, Ny, 0.0);
    Matrix Y(Nx, Ny, 0.0);

    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            X(i, j) = dx * static_cast<double>(i);
            Y(i, j) = dy * static_cast<double>(j);
        }
    }

    // Time stepping
    double dt = 2.0e-4;
    double t_end = 2.0;

    // Tolerance and logging
    double rtol = 1.0e-5;
    bool logs = true;

    // Flow field initial condition
    FlowField2D ff{
        Matrix(Nx, Ny, 0.0), // u
        Matrix(Nx, Ny, 0.0), // v
        Matrix(Nx, Ny, 0.0), // p
    };

    // REFACTOR TO HAVE 2 SEPERATE LOOPS DEPENDENT ON NX AND NY!!!!
    // Apply boundary conditions (Dirichlet)
    for (int i = 0; i < Nx; ++i) {
        ff.u(i, 0) = 0.0;
        ff.u(0, i) = 0.0;
        ff.u(Nx-1, i) = 0.0;

        ff.v(i, Ny-1) = 0.0;
        ff.v(i, 0) = 0.0;
        ff.v(0, i) = 0.0;
        ff.v(Nx-1, i) = 0.0;
        ff.u(i, Ny-1) = 1.0;
    }

    auto start = std::chrono::high_resolution_clock::now();

    solve_NavierStokes(ff, rho, nu, dt, t_end, dx, dy, Nx, Ny, rtol, logs);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Total time taken for calculations: " << duration << std::endl;
}

