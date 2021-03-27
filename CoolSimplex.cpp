//
// Created by vetinari on 14.12.20.
//

#include "CoolPoint.h"
#include "CoolConst.h"
#include "CoolSimplex.h"
#include <iostream>
#include <cmath>

void Simplex::calculate_midpoints() {
    /**
     * Calculates the coordinates of the midpoints of each face. The ith midpoint belongs to the ith face, opposite
     * of the ith vertex. It is the average of the set of vertices excluding the ith.
     */

    for (int i = 0; i < D+1; i++) { // D+1 midpoints/faces/vertices
        // Initialize current midpoint to 0
        for (int j = 0; j < D; j++) {
            midpoints[i][j] = 0;
        }

        for (int j = 0; j < D; j++) { // D *other* midpoints/faces/vertices
            int k;
            if (j < i) {
                k = j;
            } else {
                k = j+1;
            }

            for (int l = 0; l < D; l++) {   // D coordinates
                midpoints[i][l] += points[k]->coords[l];
            }
        }

        for (int j = 0; j < D; j++) {
            midpoints[i][j] /= D;
        }
    }
}


int Simplex::calculate_normals() {
    /**
     * Calculates the normals for the simplex. The ith normal is opposite the ith vertex. All normals face "outwards",
     * away from the centroid.
     *
     * return       Number of errors
     */

    int errors = 0;
    for (int i = 0; i < D+1; i++) { // Each face/point opposite
        // Build face out of all points except the one opposite
        Point * face[D];
        for (int j = 0; j < D; j++) {
            if (j < i) {
                face[j] = points[j];
            } else if (j >= i) {
                face[j] = points[j+1];
            }
        }

        // Find normal on that face
        double * norm = find_normal(&face[0]);

        double diff[D];
        double difflen2 = 0;
        for (int j = 0; j < D; j++) {
            diff[j] = centroid[j] - midpoints[i][j];
            difflen2 += diff[j] * diff[j];
        }

        double dot = 0;
        // Note: norm is normalized, diff is not
        for (int j = 0; j < D; j++) {
            dot += diff[j] * norm[j] / difflen2;
        }

        // Ensure vector points away from centroid
        if (dot > DIP_EPSILON) {
            for (int j = 0; j < D; j++) {
                norm[j] *= -1;
            }
        } else if (dot < -DIP_EPSILON) {
            // fine
        } else {
            // Count multiple errors in this one simplex as one error
            errors = 1;

#ifndef DIP_SUPPRESS_SIMPLEX_ERRORS
            std::cerr << "DIP ERROR in Simplex object " << this << " in function calculate_normals()" << std::endl;
            std::cerr << "Problem in normal calculation: Normal vector perpendicular on centroid-midpoint connection" << std::endl;
            std::cerr << "This simplex is likely degenerate, or very close to it" << std::endl;

            std::cerr << "Centroid:\t\t";
            for (int k = 0; k < D; k++) {
                std::cerr << centroid[k] << " ";
            }
            std::cerr << std::endl;

            std::cerr << "Face Midpoint:\t";
            for (int k = 0; k < D; k++) {
                std::cerr << midpoints[i][k] << " ";
            }
            std::cerr << std::endl;

            std::cerr << "Difference vector:\t" << std::endl;
            for (int j = 0; j < D; j++) {
                std::cerr << diff[j] / difflen2 << " ";
            }
            std::cerr << std::endl;

            std::cerr << "Normal vector:" << std::endl;
            for (int k = 0; k < D; k++) {
                std::cerr << norm[k] << " ";
            }
            std::cerr << std::endl;

            std::cerr << "Dot product between difference vector and normal vector: " << dot << std::endl;

            print_error_info();
#endif
        }

        // Finally, write correct normal into normals array
        for (int j = 0; j < D; j++) {
            normals[i][j] = norm[j];
        }

        delete[] norm;
    }

    return errors;
}


double * Simplex::find_normal(Point ** vertices) {
    /**
     * Finds the normal on a hyperplane defined by a set of points. Will crash if they are linearly dependent.
     *
     * The normal is given by the nullspace of the matrix of the in-plane vectors:
     * https://en.wikipedia.org/wiki/Normal_(geometry)#Hypersurfaces_in_n-dimensional_space
     *
     * The nullspace/kernel can be calculated using Gaussian elimination:
     * https://en.wikipedia.org/wiki/Kernel_(linear_algebra)#Computation_by_Gaussian_elimination
     *
     */

    double * norm = new double[D];

    // Construct matrix
    double matrix[D][D-1];      // D-1 points, D coordinates; transposed
    for (int i = 0; i < D-1; i++) {
        for (int j = 0; j < D; j++) {
            matrix[j][i] = vertices[i]->coords[j] - vertices[D-1]->coords[j];
        }
    }

    // unit matrix
    double unit[D][D];  // unit matrix
    for (int i = 0; i < D; i++) {
        for (int j = 0; j < D; j++) {
            if (i == j) {
                unit[i][j] = 1;
            } else {
                unit[i][j] = 0;
            }
        }
    }

    gauss_elimination<D, D-1, D>(matrix, unit, 0);

    int found;
    for (int i = 0; i < D; i++) {
        // If current row is full of zeros -> found norm
        int allzero = 1;
        for (int j = 0; j < D-1; j++) {
            if (matrix[i][j] != 0) {
                allzero = 0;
                break;
            }
        }

        if (allzero == 1) {
            found = 1;
            for (int j = 0; j < D; j++) {
                norm[j] = unit[i][j];
            }
            break;
        }
    }
    if (found != 1) {
        std::cerr << "Error while trying to find normal vector: Did not find zero row" << std::endl;
        std::cerr << "https://en.wikipedia.org/wiki/Kernel_(linear_algebra)#Computation_by_Gaussian_elimination" << std::endl;
        std::cerr << "Input points:" << std::endl;
        for (int i = 0; i < D; i++) {
            for (int j = 0; j < D; j++) {
                std::cerr << vertices[i]->coords[j] << " ";
            }
            std::cerr << std::endl;
        }
        std::cerr << "Initial Matrix:" << std::endl;
        for (int i = 0; i < D; i++) {
            for (int j = 0; j < D-1; j++) {
                std::cerr << vertices[j]->coords[i] - vertices[D-1]->coords[i] << "\t";
            }
            std::cerr << "|\t";
            for (int j = 0; j < D; j++) {
                if (i == j) {
                    std::cerr << 1 << "\t";
                } else {
                    std::cerr << 0 << "\t";
                }
            }
            std::cerr << std::endl;
        }
        std::cerr << "Final Matrix:" << std::endl;
        for (int i = 0; i < D; i++) {
            for (int j = 0; j < D-1; j++) {
                std::cerr << matrix[i][j] << "\t";
            }
            std::cerr << "|\t";
            for (int j = 0; j < D; j++) {
                std::cerr << unit[i][j] << "\t";
            }
            std::cerr << std::endl;
        }
        abort();
    }

    // Normalize norm
    double len = 0;
    for (int i = 0; i < D; i++) {
        len += pow(norm[i], 2);
    }
    len = sqrt(len);
    for (int i = 0; i < D; i++) {
        norm[i] /= len;
    }

    return norm;
}


template<int n>
double laplace_expansion(double (&matrix)[n][n]) {
    /**
     * Recursive function implementing Laplace expansion to calculate determinant of a matrix.
     * Slow, inefficient, simple.
     *
     * https://en.wikipedia.org/wiki/Laplace_expansion
     *
     * double * matrix      Pointer to matrix[i][j], n x n
     * int n                size of matrix
     * returns              Determinant
     */

    double det = 0;
    for (int i = 0; i < n; i++) {   // All elements of first row
        // Build submatrix without first row and ith column
        double submatrix[n-1][n-1];
        for (int j = 1; j < n; j++) {   // All rows except first
            for (int k = 0; k < n-1; k++) { // All columns except ith
                int l;
                if (k < i) {
                    l = k;
                } else {
                    l = k+1;
                }
                submatrix[j-1][k] = matrix[j][l];
            }
        }
        // The exponent is i+j. But in math indices start at one, so i becomes i+1. and j is the first row - becomes 1.
        det += pow(-1, i+1 + 1) * matrix[0][i] * laplace_expansion<n-1>(submatrix);
    }
    return det;
}


template<>
double laplace_expansion<1>(double (&matrix)[1][1]) {
    return matrix[0][0];
}



void Simplex::validate_simplex() {
    /**
     * Check that simplex is not degenerate - i.e. actualy n-d object
     * TODO doesnt work right
     */

    // Get matrix of difference vectors
    double matrix[D][D];
    for (int i = 0; i < D; i++) {
        for (int j = 0; j < D; j++) {
            matrix[j][i] = points[i]->coords[j] - points[D]->coords[j];
        }
    }

    double det = ::laplace_expansion<D>(matrix); // hacks

    if (fabs(det) < DIP_EPSILON) {
        std::cerr << "Warning in Simplex validation: Points might not be linearly independent" << std::endl;
        std::cerr << "Points:" << std::endl;
        for (int j = 0; j < D+1; j++) {
            std::cerr << "[";
            for (int k = 0; k < D; k++) {
                std::cerr << points[j]->coords[k] << ",\t";
            }
            std::cerr << "], " << std::endl;
        }
        std::cerr << "Matrix of difference vectors: " << std::endl;
        for (int j = 0; j < D; j++) {
            for (int k = 0; k < D; k++) {
                std::cerr << matrix[j][k] << "\t";
            }
            std::cerr << std::endl;
        }
        std::cerr << "Determinant: " << det << std::endl;
    }
}


int Simplex::validate_normals() {
    /**
     * Check that all normals are
     * a) Normalized
     * b) Normal
     * c) Pointing away from the centroid
     * d) Pointing away from each other
     *
     * return       0 if no errors, 1 if errors
     */

    int error = 0;

    // a) Validate normalisation
    double normfailed[D+1] = { 0 }; // Array of zeros
    int nf_flag = 0;
    for (int i = 0; i < D+1; i++) {
        double len2 = 0;
        for (int j = 0; j < D; j++) {
            len2 += pow(normals[i][j], 2);
        }

        if (fabs(len2 - 1) > DIP_EPSILON) {
            nf_flag = 1;
            normfailed[i] = sqrt(len2);
        }
    }
    if (nf_flag) {
#ifndef DIP_SUPPRESS_SIMPLEX_ERRORS
        std::cerr << "DIP ERROR in Simplex object " << this << " in function validate_normals()" << std::endl;
        std::cerr << "Error in normal validation: Normalisation of normals " << std::endl;
#endif
        error = 1;
    }

    // b) Normal on hyperplanes
    for (int i = 0; i < D+1; i++) { // all faces/normals
        Point * facepoints[D];
        for (int j = 0; j < D; j++) {   // All points in face/All points in simplex except point i
            int k;
            if (j < i) {
                k = j;
            } else {
                k = j + 1;
            }

            facepoints[j] = points[k];
        }

        for (int j = 1; j < D; j++) {
            double diff[D]; // In the hyperplane
            double len = 0;

            for (int k = 0; k < D; k++) {
                diff[k] = facepoints[0]->coords[k] - facepoints[j]->coords[k];
                len += pow(diff[k], 2);
            }
            len = sqrt(len);
            double dot = 0;
            for (int k = 0; k < D; k++) {
                dot += (diff[k] / len) * normals[i][k];
            }

            if (fabs(dot) > DIP_EPSILON) {
#ifndef DIP_SUPPRESS_SIMPLEX_ERRORS
                std::cerr << "DIP ERROR in Simplex object " << this << " in function validate_normals()" << std::endl;
                std::cerr << "Normal vector not normal: " << i << std::endl;
                std::cerr << "Normal vector: " << std::endl;
                for (int k = 0; k < D; k++) {
                    std::cerr << normals[i][k] << " ";
                }
                std::cerr << std::endl;
                std::cerr << "In plane vector: " << std::endl;
                for (int k = 0; k < D; k++) {
                    std::cerr << diff[k] / len << " ";
                }
                std::cerr << std::endl;
                std::cerr << "Dot product: " << dot << std::endl;
                std::cerr << "Points used for in plane vector: " << 0 << " " << j << std::endl;
                std::cerr << "Plane points:" << std::endl;

                for (int k = 0; k < D; k++) {
                    std::cerr << k << ":\t";
                    for (int l = 0; l < D; l++) {
                        std::cerr << facepoints[k]->coords[l] << " ";
                    }
                    std::cerr << std::endl;
                }
#endif
                error = 1;
            }

        }
    }

    // c)
    // Pointing away from centroid - use scalar product between normal and difference vector pointing from midpoint
    // to centroid
    for (int i = 0; i < D+1; i++) {
        double diff[D];
        double len = 0;
        for (int j = 0; j < D; j++) {
            diff[j] = centroid[j] - midpoints[i][j];
            len += pow(diff[j], 2);
        }
        len = sqrt(len);

        double dot = 0;
        for (int j = 0; j < D; j++) {
            dot += normals[i][j] * (diff[j] / len);
        }
        if (dot > DIP_EPSILON) {
#ifndef DIP_SUPPRESS_SIMPLEX_ERRORS
            std::cerr << "DIP ERROR in Simplex object " << this << " in function validate_normals()" << std::endl;
            std::cerr << "Normal vector not pointing away from centroid: " << i << std::endl;
            for (int j = 0; j < D; j++) {
                std::cerr << normals[i][j] << " ";
            }
            std::cerr << std::endl;
            std::cerr << "Midpoint vector:" << std::endl;
            for (int j = 0; j < D; j++) {
                std::cerr << midpoints[i][j] << " ";
            }
            std::cerr << std::endl;
            std::cerr << "Difference vector: " << std::endl;
            for (int j = 0; j < D; j++) {
                std::cerr << diff[j] << " ";
            }
            std::cerr << std::endl;
            std::cerr << "Difference vector normalized: " << std::endl;
            for (int j = 0; j < D; j++) {
                std::cerr << diff[j] / len << " ";
            }
            std::cerr << std::endl;
            std::cerr << "Dot product: " << dot << std::endl;
#endif
            error = 1;
        }
    }

    // d)
    // The largest angle is the smallest dot product
    // This dot product cannot be larger than 0 - the two normals with the largest angle cannot point in "roughly the
    // same direction"
    // In fact, the "smallest possible largest angle" is given for a regular simplex (equilateral triangle, tetrahedron, ...)
    // This angle is given by arccos(-1/D)
    // Use this criterium for a sanity check.
    double smallest_dot = 1;
    int smallest_dot_idx1, smallest_dot_idx2;
    for (int i = 0; i < D+1; i++) {
        for (int j = 0; j < D+1; j++) {
            double dot = 0;
            for (int k = 0; k < D; k++) {
                dot += normals[i][k] * normals[j][k];
            }
//            std::cout << i << " " << j << " " << dot << std::endl;
//            std::cout << "i: " << i << " j: " << j << " dot: " << dot << std::endl;
            if (dot < smallest_dot) {
                smallest_dot = dot;
                smallest_dot_idx1 = i;
                smallest_dot_idx2 = j;
            }
        }
    }
//    std::cout << std::endl;

    if (smallest_dot > -1./D) {
#ifndef DIP_SUPPRESS_SIMPLEX_ERRORS
        std::cerr << "DIP ERROR in Simplex object " << this << " in function validate_normals()" << std::endl;
        std::cerr << "Error in normal validation: Scalar product of Normals "
                  << smallest_dot_idx1 << " and " << smallest_dot_idx2 << ": "
                  << smallest_dot << " > -1/" << D << " = " << -1./D << std::endl;
#endif
        error = 1;
    }

#ifndef DIP_SUPPRESS_SIMPLEX_ERRORS
    if (error == 1) {
        print_error_info();
    }
#endif

    return error;
}


void Simplex::construct_T_inv() {
    /**
     * Construct the matrix required to calculate the barycentric coordinates relative to this simplex
     * https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Barycentric_coordinates_on_tetrahedra
     */
    calculate_T();
    invert_T();
}


void Simplex::calculate_T() {
    for (int i = 0; i < D; i++) {   // for each point (except the last)
        for (int j = 0; j < D; j++) {   // for each coordinate
            // matrix[j][i], not matrix[i][j] - coordinates go down, points go right
            T_inv[j][i] = points[i]->coords[j] - points[D]->coords[j];
        }
    }
}


void Simplex::invert_T() {
    double unit[D][D];  // unit matrix
    for (int i = 0; i < D; i++) {
        for (int j = 0; j < D; j++) {
            if (i == j) {
                unit[i][j] = 1;
            } else {
                unit[i][j] = 0;
            }
        }
    }

    gauss_elimination<D, D, D>(T_inv, unit, 1);
//    for (int i = 0; i < D; i++) {   // for each point (except the last)
//        for (int j = 0; j < D; j++) {   // for each coordinate
//             matrix[j][i], not matrix[i][j] - coordinates go down, points go right
//            T[i][j] = T_inv[i][j];
//        }
//    }

    for (int i = 0; i < D; i++) {
        for (int j = 0; j < D; j++) {
            T_inv[i][j] = unit[i][j];
        }
    }
}



template<int M, int N1, int N2>
void Simplex::gauss_elimination(double (&A)[M][N1], double (&B)[M][N2], int apply_permutation) {
    /**
     * Ax = B
     * Performs row wise Gauss Jordan elimination on A with given dimensions. Transformations are also applied to B.
     * Uses partial pivoting. Adapted from Numerical Recipes chapter 2.1
     *
     * Works in place on A and B.
     *
     * TODO Implicit pivoting
     *
     * double A[M][N1]          The matrix to solve
     * double B[M][N2]          The other matrix to apply the transformations to; can also be a solution vector etc.
     * int apply_permutation    Whether to apply the row permutations to the matrices in the end; important for some
     *                          applications (inverse calculation), irrelevant for others (calculating normals).
     *                          0 for no, 1 for yes
     */

    // Permutation of ROW InDeXes
    int row_idx[D];
    for (int i = 0; i < M; i++) {
        row_idx[i] = i;
    }

    for (int i = 0; i < M; i++) {
        // find pivot
        int pivot_row_idx = i;
        double pivot = A[row_idx[i]][i];
        for (int j = i; j < M; j++) {
            if (fabs(A[row_idx[j]][i]) > fabs(pivot)) {
                pivot = A[row_idx[j]][i];
                pivot_row_idx = j;
            }
        }
        if (fabs(pivot) < DIP_EPSILON) {
            continue;
        }

        // pivoting
        int temp = row_idx[i];
        row_idx[i] = row_idx[pivot_row_idx];
        row_idx[pivot_row_idx] = temp;

        // Normalization of current row
        double norm_factor = A[row_idx[i]][i];
        for (int j = 0; j < N1; j++) {
            A[row_idx[i]][j] /= norm_factor;
        }
        for (int j = 0; j < N2; j++) {
            B[row_idx[i]][j] /= norm_factor;
        }

        // Subtracting the current row from all other rows
        for (int j = 0; j < M; j++) {
            if (j != i) {
                double factor = A[row_idx[j]][i];
                for (int k = 0; k < N1; k++) {
                    A[row_idx[j]][k] -= factor * A[row_idx[i]][k];
                }
                for (int k = 0; k < N2; k++) {
                    B[row_idx[j]][k] -= factor * B[row_idx[i]][k];
                }
            }
        }

    }

    if (apply_permutation) {
        for (int i = 0; i < M; i++) {
            // Switch elements
            for (int j = 0; j < N1; j++) {
                double temp = A[row_idx[i]][j];
                A[row_idx[i]][j] = A[i][j];
                A[i][j] = temp;
            }
            for (int j = 0; j < N2; j++) {
                double temp = B[row_idx[i]][j];
                B[row_idx[i]][j] = B[i][j];
                B[i][j] = temp;
            }

            // Switch index. These lines cost eight hours of my life.
            for (int j = 0; j < M; j++) {
                if (row_idx[j] == i) {
                    row_idx[j] = row_idx[i];
                    break;
                }
            }
            row_idx[i] = i;
        }
    }
}


double * Simplex::convert_to_bary(const double * coords) {
    /**
     * Convert the given coordinates to the barycentric coordinate system of the simplex.
     *
     * https://math.stackexchange.com/a/1226825
     */
    // 1. Obtain p - v_n+1
    double vec[D];
//    std::cout << "T: " << std::endl;
//    for (int i = 0; i < D; i++) {
//        for (int j = 0; j < D; j++) {
//            std::cout << T[i][j] << " ";
//        }
//        std::cout << std::endl;
//    }
//    std::cout << "coord: " << coords[0] << " " << coords[1] << " " << coords[2] << std::endl;
//    std::cout << "p: " << points[D]->coords[0] << " " << points[D]->coords[1] << " " << points[D]->coords[2] << std::endl;
//    std::cout << "coord - p: ";
    for (int i = 0; i < D; i++) {
        vec[i] = coords[i] - points[D]->coords[i];
//        std::cout << vec[i] << " ";
    }
//    std::cout << std::endl;

    // 2. Multiply T_inv * (p - v_n+1) to get lambda
    double * bary = new double[D+1];
    for (int i = 0; i < D; i++) {
        bary[i] = 0;
        for (int j = 0; j < D; j++) {
            bary[i] += T_inv[i][j] * vec[j];
//            std::cout << T_inv[i][j] << " * " << vec[j] << " + ";
        }
//        std::cout << " = " << bary[i] << std::endl;
    }

    // dependent coordinate lambda_n+1
    bary[D] = 1;
    for (int i = 0; i < D; i++) {
        bary[D] -= bary[i];
    }

    return bary;
}


int Simplex::check_bary(const double* bary) {
    /**
     * Check if the given barycentric coordinates belong to a point inside or outside of the simplex
     * Returns 1 if inside, 0 otherwise.
     * Edge cases are considered to be outside.
     */
    double bsum = 0;
    for (int i = 0; i < D+1; i++) {
        if (bary[i] < 0) {
            // Consider points on edges to not be contained, in order to avoid degeneracies!
            return 0;
        }
        if (isnanf(bary[i])) {
            std::cerr << "DIP ERROR in Simplex object " << this << " in function check_bary()" << std::endl;
            std::cerr << "Error: Barycentric coordinate " << i << " is nan" << std::endl;
            std::cerr << "Something has gone SERIOUSLY wrong!" << std::endl;
            print_error_info();
            abort();
        }
        bsum += bary[i];
    }

    if (bsum > 1 + DIP_EPSILON) {
        return 0;
    }

    return 1;
}


void Simplex::print_error_info() {
    /**
     * Print general information about this simplex to std::cerr.
     */
    std::cerr << "======================" << std::endl;
    std::cerr << "Error Info of Simplex " << index << " at address " << this << ":" << std::endl;
    std::cerr << "Centroid:" << std::endl;
    for (int i = 0; i < D; i++) {
        std::cerr << centroid[i] << " ";
    }
    std::cerr << std::endl;

    std::cerr << "Indices of neighbor simplices:" << std::endl;
    for (int i = 0; i < D+1; i++) {
        std::cerr << neighbour_indices[i] << " ";
    }
    std::cerr << std::endl;

    std::cerr << "Points: " << std::endl;
    for (int i = 0; i < D+1; i++) {
        for (int j = 0; j < D; j++) {
            std::cerr << points[i]->coords[j] << " ";
        }
        std::cerr << points[i]->value << "\t(" << points[i] << ")" << std::endl;
    }

    std::cerr << "Midpoints of faces: " << std::endl;
    for (int i = 0; i < D+1; i++) {
        for (int j = 0; j < D; j++) {
            std::cerr << midpoints[i][j] << " ";
        }
        std::cerr << std::endl;
    }

    std::cerr << "Normal vectors on faces:" << std::endl;
    for (int i = 0; i < D+1; i++) {
        for (int j = 0; j < D; j++) {
            std::cerr << normals[i][j] << " ";
        }
        std::cerr << std::endl;
    }

    std::cerr << "T_inv:" << std::endl;
    for (int i = 0; i < D; i++) {
        for (int j = 0; j < D; j++) {
            std::cerr << T_inv[i][j] << " ";
        }
        std::cerr << std::endl;
    }
    std::cerr << "======================" << std::endl << std::endl;
}