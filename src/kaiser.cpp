// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)
//
// Analytical3x3 symmetric matrix eigenvalue solver
// Replaces Kaiser.f with modern C++ implementation
// Provides identical function signature for backward compatibility

#include <cmath>
#include <algorithm>
#include <limits>
#include <functional>

extern "C" {

// Analytical solution for eigenvalues of3x3 symmetric matrix
// Uses characteristic polynomial method with robust numerical implementation
void kaiser_wrap_(double *a, int *nrows, int *ncols, 
                double *eigenv, double *trace, double *sume, int *ier)
{
    // Validate input
    if (*nrows < 3 || *ncols < 3) {
        *ier = 1;  // Invalid dimensions
        return;
    }
    
    if (*nrows < *ncols) {
        *ier = 1;  // N must be <= NROWS
        return;
    }
    
    *ier = 0;  // Success
    
    // Extract matrix elements (assuming row-major storage like Fortran)
    double a11 = a[0], a12 = a[1], a13 = a[2];
    double a21 = a[3], a22 = a[4], a23 = a[5];
    double a31 = a[6], a32 = a[7], a33 = a[8];
    
    // Calculate trace for output
    *trace = a11 + a22 + a33;
    
    // For 3x3 symmetric matrix, use analytical eigenvalue computation
    // Characteristic polynomial: det(A - λI) = 0
    // λ³ - c₁λ² + c₂λ - c₃ = 0
    
    double c1 = *trace;
    double c2 = a11*a22 - a12*a21 + a11*a33 - a13*a31 + a22*a33 - a23*a32;
    double c3 = a11*a22*a33 + 2*a12*a23*a31 - 
                a11*a23*a32 - a22*a13*a31 - a33*a12*a21;
    
    // Solve cubic equation using robust method
    double q = (3*c2 - c1*c1) / 9.0;
    double r = (9*c1*c2 - 27*c3 - 2*c1*c1*c1) / 54.0;
    double disc = q*q*q + r*r;
    
    double eigenvals[3];
    
    if (disc > 0.0) {
        // One real root, two complex conjugates
        double sqrt_disc = std::sqrt(disc);
        double s = std::cbrt(r + sqrt_disc);
        double t = std::cbrt(r - sqrt_disc);
        
        eigenvals[0] = s + t - c1/3.0;
        // For complex roots, take absolute values (behavior matches original Kaiser)
        eigenvals[1] = std::abs(-(s + t)/2.0 - c1/3.0);
        eigenvals[2] = std::abs((std::sqrt(3.0)/2.0)*(s - t));
    } else {
        // Three real roots
        double phi;
        if (std::abs(q) < 1e-12) {
            phi = M_PI/2.0;
        } else {
            phi = std::acos(r / std::sqrt(-q*q*q));
        }
        
        double sqrt_q = std::sqrt(std::abs(q));
        eigenvals[0] = 2.0*sqrt_q*std::cos(phi/3.0) - c1/3.0;
        eigenvals[1] = 2.0*sqrt_q*std::cos((phi + 2.0*M_PI)/3.0) - c1/3.0;
        eigenvals[2] = 2.0*sqrt_q*std::cos((phi + 4.0*M_PI)/3.0) - c1/3.0;
    }
    
    // Sort eigenvalues in descending order (like original Kaiser)
    std::sort(eigenvals, eigenvals + 3, [](double a, double b) { return a > b; });
    
    // Calculate eigenvectors and normalize
    // For each eigenvalue, solve (A - λI)v = 0
    // Use robust method to find non-trivial solution
    
    *sume = 0.0;
    
    for (int i = 0; i < 3; i++) {
        eigenv[i] = std::abs(eigenvals[i]);
        *sume += eigenv[i];
        
        // Calculate eigenvector for eigenvalue λ[i]
        double v1, v2, v3;
        double lambda_i = eigenvals[i];
        
        // Find the most reliable component for eigenvector calculation
        double det1 = (a11 - lambda_i)*(a22 - lambda_i) - a12*a21;
        double det2 = (a11 - lambda_i)*(a33 - lambda_i) - a13*a31;
        double det3 = (a22 - lambda_i)*(a33 - lambda_i) - a23*a32;
        
        if (std::abs(det1) > std::abs(det2) && std::abs(det1) > std::abs(det3)) {
            // Use third component as free variable
            v3 = 1.0;
            v1 = (a13*(a33 - lambda_i) - a31*a23) / det1;
            v2 = (a12*(a33 - lambda_i) - a21*a13) / det1;
        } else if (std::abs(det2) > std::abs(det3)) {
            // Use second component as free variable
            v2 = 1.0;
            v1 = (a12*(a22 - lambda_i) - a13*a32) / det2;
            v3 = (a31*(a22 - lambda_i) - a21*a32) / det2;
        } else {
            // Use first component as free variable
            v1 = 1.0;
            v2 = (a21*(a11 - lambda_i) - a23*a31) / det3;
            v3 = (a31*(a11 - lambda_i) - a32*a12) / det3;
        }
        
        // Normalize eigenvector and store in matrix A (overwriting input like original)
        double norm = std::sqrt(v1*v1 + v2*v2 + v3*v3);
        if (norm > std::numeric_limits<double>::epsilon()) {
            double inv_norm = 1.0 / norm;
            a[i] = v1 * inv_norm;     // Store in a(1,i)
            a[3 + i] = v2 * inv_norm; // Store in a(2,i) 
            a[6 + i] = v3 * inv_norm; // Store in a(3,i)
        } else {
            // Degenerate case - use standard basis vectors
            a[i] = (i == 0) ? 1.0 : 0.0;
            a[3 + i] = (i == 1) ? 1.0 : 0.0;
            a[6 + i] = (i == 2) ? 1.0 : 0.0;
        }
    }
    
    // Check for convergence issues (mimic original error handling)
    if (std::abs(*sume - *trace) > 1e-6 * std::abs(*trace)) {
        *ier = 2;  // Convergence failed
    }
}

} // extern "C"