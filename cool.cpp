//
// Created by vetinari on 15.09.20.
//

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cfloat>
#include <cmath>
#include <cassert>
#include <vector>
#include <map>
#include <unordered_set>
#include <algorithm>
#include <bitset>

// TODO ?
#define EPSILON 0.001

// https://stackoverflow.com/a/4609795
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


template<int D> class Point;
template<int D> class Simplex;
template<int N, int D, int S> class Cool;

template<int D>
class Point {
    /**
     * Class representing a single data point. Mostly used to find simplices it is part of.
     */
    template<int, int, int> friend class Cool;
    template<int> friend class Simplex;
private:
    double value;

public:
    double coords[D];
};


template<int D>
class Simplex {
    /**
     * Class representing an M-d simplex in the triangulation. Has M+1 vertices.
     *
     * int M        Number of dimensions
     */
    template<int, int, int> friend class Cool;
private:
    // TODO I would like to make these const, but I don't think I can, since the value is determined at runtime
    Point<D> * points[D+1];             // D+1 points; Array of pointers to Point<D>
    double T_inv[D][D];                 // T: https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Barycentric_coordinates_on_tetrahedra
    int neighbour_indices[D+1];         // One neighbour opposite every point
    Simplex * neighbour_pointers[D+1];
    double centroid[D];
    double sbtree_radius_sq;
    double midpoints[D+1][D];           // Midpoints of the faces; D+1 faces, D coordinates
    double normals[D+1][D];             // Outward pointing normals of the faces; D+1 faces, D coordinates

    void invert_T();
    double * convert_to_bary(const double *);
    int check_bary(const double *);
    void calculate_normals();
    void calculate_midpoints();

    template<size_t M, size_t N1, size_t N2> void gauss_elimination(double (&)[M][N1], double (&)[M][N2], int);
public:
    double * find_normal(Point<D> ** );
    void validate_simplex();
    void validate_normals();
    Simplex * lchild;
    Simplex * rchild;

    Simplex() {
        sbtree_radius_sq = 0;
    };

    void calculate_centroid() {
        /**
         * Calculate centroid (=avg) from points.
         */
        for (int i = 0; i < D; i++) {   // coordinates
            centroid[i] = 0;
            for (int j = 0; j < D+1; j++) {     // points.csv
                centroid[i] += points[j]->coords[i];
            }
            centroid[i] /= (D+1);
        }
    }
    template<size_t n> double laplace_expansion(double (&)[n][n]);
};

template<int N, int D, int S>
class Cool {
    /**
     * int N        Number of points.csv
     * int D        Dimensionality of points.csv
     * int S        Number of Simplices
     */
private:
    inline static Point<D> points[N];
    inline static Simplex<D> simplices[S];
    Simplex<D> * sbtree;         // Points to the root of the simplex ball tree

    Simplex<D> * construct_simplex_btree_recursive(Simplex<D> **, int);
    Simplex<D> * find_nearest_neighbour_sbtree(Simplex<D> *, const double *, Simplex<D> *, double);

    int sbtree_flips, sbtree_interpolate_calls;
public:
    Cool() {
        sbtree_flips = 0;
        sbtree_interpolate_calls = 0;
        avg_sbtree_flips = 0;
        for (int i = 0; i < D; i++) {
            mins[i] =      DBL_MAX;
            maxs[i] = -1 * DBL_MAX;
        }
    };
    double avg_sbtree_flips;

    // Minimum and maximum values in each dimension
    double mins[D];
    double maxs[D];

    int read_files(std::string, std::string, std::string);
    void save_sbtree(std::string);
    int construct_simplex_btree();
    double interpolate_sbtree(double *);
};


template<int D>
void Simplex<D>::calculate_normals() {
    /**
     * Calculates the normals for the simplex. The ith normal is opposite the ith vertex. All normals face "outwards",
     * away from the centroid.
     */

//    std::cout << "Calculating normals for:" << std::endl;
//    for (int i = 0; i < D+1; i++) {
//        for (int j = 0; j < D; j++) {
//            std::cout << points[i]->coords[j] << " ";
//        }
//        std::cout << std::endl;
//    }
//    std::cout << "Centroid:" << std::endl;
//    for (int i = 0; i < D; i++) {
//        std::cout << centroid[i] << " ";
//    }
//    std::cout << std::endl;
    for (int i = 0; i < D+1; i++) { // Each face/point opposite
        std::cout << "\nNormal " << i << std::endl;
        Point<D> * face[D];

        // All points except the one opposite
        for (int j = 0; j < D; j++) {
            if (j < i) {
                face[j] = points[j];
            } else if (j >= i) {
                face[j] = points[j+1];
            }
        }
//        std::cout << std::endl << "Face " << i << std::endl;
//        for (int j = 0; j < D; j++) {
//            for (int k = 0; k <D; k++) {
//                std::cout << face[j]->coords[k] << " ";
//            }
//            std::cout << std::endl;
//        }

        double * norm = find_normal(&face[0]);

        std::cout << "Found normal: " << std::endl;
        for (int j = 0; j < D; j++) {
            std::cout << norm[j] << " ";
        }
        std::cout << std::endl;

        // Make sure norm points away from centroid
        // Go +1/-1 in direction of the normal from the centroid.
        // If +1 is closer to the midpoint, normal points in the right direction. Else, flip
        double step_pos[D];
        double step_neg[D];
        double dist2_pos = 0;
        double dist2_neg = 0;

        for (int j = 0; j < D; j++) {
            step_pos[j] = centroid[j] + norm[j];
            step_neg[j] = centroid[j] - norm[j];
        }

        std::cout << step_pos[0] << " " << step_pos[1] << std::endl;
        std::cout << step_neg[0] << " " << step_neg[1] << std::endl;
        std::cout << midpoints[i][0] << " " << midpoints[i][1] << std::endl;

        for (int j = 0; j < D; j++) {
            dist2_pos += pow(step_pos[j] - midpoints[i][j], 2);
            dist2_neg += pow(step_neg[j] - midpoints[i][j], 2);
        }

        std::cout << "pos dist: " << sqrt(dist2_pos) << std::endl;
        std::cout << "neg dist: " << sqrt(dist2_neg) << std::endl;

        if (dist2_neg < dist2_pos) {
            std::cout << "Flip" << std::endl;
            for (int j = 0; j < D; j++) {
                norm[j] *= -1;
            }
        }

        // Finally, write correct normal into normals array
        for (int j = 0; j < D; j++) {
            normals[i][j] = norm[j];
        }

        delete[] norm;
    }

    std::cout << std::endl;
}


template<int D>
double * Simplex<D>::find_normal(Point<D> ** vertices) {
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
//            std::cout << "Zero row: " << i << std::endl;
            found = 1;
            for (int j = 0; j < D; j++) {
                norm[j] = unit[i][j];
//                std::cout << norm[j] << " ";
            }
//            std::cout << std::endl;
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
//    std::cout << "Len: " << len << std::endl;
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

/*
template<>
double laplace_expansion<1>(double (&matrix)[1][1]) {
    return matrix[0][0];
}
 */


template<int D>
void Simplex<D>::validate_simplex() {
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

    if (fabs(det) < EPSILON) {
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


template<int D>
void Simplex<D>::validate_normals() {
    /**
     * Check that all normals are
     * a) Normalized
     * b) Normal
     * c) Pointing away from the centroid
     * d) Pointing away from each other
     */

    // a) Validate normalisation
    double normfailed[D+1] = { 0 }; // Array of zeros
    int nf_flag = 0;
    for (int i = 0; i < D+1; i++) {
        double len2 = 0;
        for (int j = 0; j < D; j++) {
            len2 += pow(normals[i][j], 2);
        }

        if (fabs(len2 - 1) > EPSILON) {
            nf_flag = 1;
            normfailed[i] = sqrt(len2);
        }
    }
    if (nf_flag) {
        std::cerr << "Error in normal validation: Normalisation of normals " << std::endl;

        std::cerr << "Simplex points:" << std::endl;
        for (int j = 0; j < D+1; j++) {
            std::cerr << j << ":   ";
            for (int k = 0; k < D; k++) {
                std::cerr << points[j]->coords[k] << " ";
            }
            std::cerr << std::endl;
        }
        std::cerr << "Centroid: " << std::endl;
        for (int j = 0; j < D; j++) {
            std::cerr << centroid[j] << " ";
        }
        std::cerr << std::endl;
        std::cerr << "Normal vectors:" << std::endl;
        for (int j = 0; j < D+1; j++) {
            if (normfailed[j] != 0) {
                std::cerr << ">" << j << ":  ";
            } else {
                std::cerr << j << ":   ";
            }
            for (int k = 0; k < D; k++) {
                std::cerr << normals[j][k] << " ";
            }
            if (normfailed[j] != 0) {
                std::cerr << "\t--> " << normfailed[j];
            }
            std::cerr << std::endl;
        }
        abort();
    }

    // b) Normal on hyperplanes
    for (int i = 0; i < D+1; i++) { // all faces/normals
        Point<D> * facepoints[D];
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

            if (fabs(dot) > EPSILON) {
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
                abort();
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
        if (dot > 0) {
            std::cerr << "Normal vector not pointing away from centroid: " << i << std::endl;
            std::cerr << "Normal vector:" << std::endl;
            for (int j = 0; j < D; j++) {
                std::cerr << normals[i][j] << " ";
            }
            std::cerr << std::endl;
            std::cerr << "Centroid: " << std::endl;
            for (int j = 0; j < D; j++) {
                std::cerr << centroid[j] << " ";
            }
            std::cerr << "Midpoint vector:" << std::endl;
            for (int j = 0; j < D; j++) {
                std::cerr << midpoints[i][j] << " ";
            }
            std::cerr << "Difference vector: " << std::endl;
            for (int j = 0; j < D; j++) {
                std::cerr << diff[j] << " ";
            }
            std::cerr << "Dot product: " << dot << std::endl;
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
//            std::cout << "i: " << i << " j: " << j << " dot: " << dot << std::endl;
            if (dot < smallest_dot) {
                smallest_dot = dot;
                smallest_dot_idx1 = i;
                smallest_dot_idx2 = j;
            }
        }
    }

    if (smallest_dot > -1./D) {
        std::cerr << "Error in normal validation: Scalar product of Normals "
                  << smallest_dot_idx1 << " and " << smallest_dot_idx2 << ": "
                  << smallest_dot << " < 1/" << D << " = " << 1./D << std::endl;

        std::cerr << "Simplex points:" << std::endl;
        for (int j = 0; j < D+1; j++) {
            std::cerr << j << ":   ";
            for (int k = 0; k < D; k++) {
                std::cerr << points[j]->coords[k] << " ";
            }
            std::cerr << std::endl;
        }
        std::cerr << "Centroid: " << std::endl;
        for (int j = 0; j < D; j++) {
            std::cerr << centroid[j] << " ";
        }
        std::cerr << std::endl;
        std::cerr << "Midpoints:" << std::endl;
        for (int j = 0; j < D+1; j++) {
            std::cerr << j << ":   ";
            for (int k = 0; k < D; k++) {
                std::cerr << midpoints[j][k] << " ";
            }
            std::cerr << std::endl;
        }
        std::cerr << "Normal vectors:" << std::endl;
        for (int j = 0; j < D+1; j++) {
            if (j == smallest_dot_idx1 || j == smallest_dot_idx2) {
                std::cerr << ">" << j << ":  ";
            } else {
                std::cerr << j << ":   ";
            }
            for (int k = 0; k < D; k++) {
                std::cerr << normals[j][k] << " ";
            }
            std::cerr << std::endl;
        }

        abort();
    }
}


template<int D>
void Simplex<D>::calculate_midpoints() {
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


template<int N, int D, int S>
int Cool<N, D, S>::construct_simplex_btree() {
    /**
     * This function serves as a public "adapter" to the actual ball tree construction function,
     * construct_point_btree_recursive().
     *
     * Returns 0 on success.
     */
    // Construct array of pointers
    Simplex<D> * simps[S];
    for (int i = 0; i < S; i++) {
        simps[i] = &(simplices[i]);
    }

    sbtree = construct_simplex_btree_recursive(simps, S);

    return 0;
}


template< int N, int D, int S>
Simplex<D> * Cool<N, D, S>::construct_simplex_btree_recursive(Simplex<D> ** simps, int n) {
    /**
     * Simplex<D> simps     Pointer to array of pointers to simplices to organize in tree
     * int n                Number of elements in array
     */

    // Input validation + what if theres only one element left?
    if (n < 0) {
        std::cerr << "Illegal argument in construct_btree: n = " << n << std::endl;
    } else if (n == 1) {
        simps[0]->lchild = NULL;
        simps[0]->rchild = NULL;
        simps[0]->sbtree_radius_sq = 0;

        return *simps;
    }

    // More than one element
    // 1. Find dimension of greatest spread
    // Working with the centroids of the triangles
    double largest_spread = -1;
    int lspread_dim = -1;
    double avg = 0;             // Traditionally Ball trees work with the median, but that requires sorting TODO

    for (int i = 0; i < D; i++) {   // every dimension
        double dim_min =    DBL_MAX;
        double dim_max = -1*DBL_MAX;
        double dim_spread = 0;
        for (int j = 0; j < n; j++) {   // every simplex
            double val = simps[j]->centroid[i];
            if (val < dim_min) {
                dim_min = val;
            }
            if (val > dim_max) {
                dim_max = val;
            }
        }
//        std::cout << "Dim: " << i << " min: " << dim_min << " max: " << dim_max << std::endl;

        assert(dim_min <= dim_max);
        dim_spread = fabs(dim_max - dim_min);
        if (dim_spread > largest_spread) {
            largest_spread = dim_spread;
            lspread_dim    = i;
            avg            = (dim_max + dim_min) / 2;
        }
    }
//    std::cout << "spread dim: " << lspread_dim << " avg: " << avg << std::endl;

    assert(largest_spread != -1 && lspread_dim != -1);

    // 2. Find pivot - simplex with centroid closest to avg          and
    // 3. Group points into sets to the left and right of average (L, R)
    double min_dist = DBL_MAX;
    Simplex<D> * pivot_addr = NULL;

    Simplex<D> ** L = new Simplex<D> * [n];
    Simplex<D> ** R = new Simplex<D> * [n];
    int lcount = 0;
    int rcount = 0;

    for (int i = 0; i < n; i++) {
        double val = simps[i]->centroid[lspread_dim];
        Simplex<D> * addr = simps[i];
        double dist = fabs(avg - val);

        if (dist < min_dist) {
            // Add old pivot to L or R
            if (pivot_addr) {
                if (pivot_addr->centroid[lspread_dim] <= avg) {
                    L[lcount] = pivot_addr;
                    lcount++;
                } else {
                    R[rcount] = pivot_addr;
                    rcount++;
                }
            }

            // Save new pivot
            min_dist = dist;
            pivot_addr = simps[i];

        } else {
            if (val <= avg) {
                L[lcount] = addr;
                lcount++;
            } else {
                R[rcount] = addr;
                rcount++;
            }
        }
    }
    assert(pivot_addr != NULL);
    assert(rcount != 0 || lcount != 0);
    assert(min_dist < DBL_MAX);
    for (int i = 0; i < lcount; i++) {
        assert(L[i] != pivot_addr);
    }
    for (int i = 0; i < rcount; i++) {
        assert(R[i] != pivot_addr);
    }

    // 4. Recurse on L and R
    if (lcount > 0) {
        pivot_addr->lchild = construct_simplex_btree_recursive(L, lcount);
        assert(pivot_addr->lchild != pivot_addr);
    } else {
        pivot_addr->lchild = NULL;
    }
    if (rcount > 0) {
        pivot_addr->rchild = construct_simplex_btree_recursive(R, rcount);
    } else {
        pivot_addr->rchild = NULL;
    }

    // This should only happen if there is only one element left, and in that case this line shouldnt be reached
    assert(pivot_addr->lchild != NULL || pivot_addr->rchild != NULL);

    // 5. Determine ball radius
    pivot_addr->sbtree_radius_sq = 0;
    for (int i = 0; i < n; i++) {
        double dist = 0;
        for (int j = 0; j < D; j++) {
            dist += pow(simps[i]->centroid[j] - pivot_addr->centroid[j], 2);
        }

        if (dist > pivot_addr->sbtree_radius_sq) {
            pivot_addr->sbtree_radius_sq = dist;
        }
    }

    delete[] L;
    delete[] R;

    return pivot_addr;
}


template<int N, int D, int S>
void Cool<N, D, S>::save_sbtree(std::string filename) {
    /**
     * Saves the Ball tree. Format: [Index of simplex] [Index of left child] [Index of right child] [pbtree_radius_sq]
     *
     * Grows with O(N^2) ... TODO
     */
    std::ofstream file;
    file.open(filename);
    for (int i = 0; i < S; i++) {
        file << i << " ";
        int lchild = -1, rchild = -1;
        for (int j = 0; j < S; j++) {
            if (simplices[i].lchild == &simplices[j]) {
                lchild = j;
            }
            if (simplices[i].rchild == &simplices[j]) {
                rchild = j;
            }
        }
        file << lchild << " " << rchild << " " << simplices[i].sbtree_radius_sq << std::endl;
    }
    file.close();
}


template<int N, int D, int S>
Simplex<D> * Cool<N, D, S>::find_nearest_neighbour_sbtree(Simplex<D> * root, const double * target, Simplex<D> * best, double min_dist2) {
    /**
     * Recursive function to find the nearest neighbor of point target in tree root.
     * Adapted from https://en.wikipedia.org/wiki/Ball_tree#Pseudocode_2
     *
     * TODO I think this algorithm can be optimized to reduce the number of distance calculations
     *  -> talk to Klaus
     *
     * root             Root of subtree to search in
     * target           Coordinates of input point
     * best             Simplex with closest centroid found so far
     * min_dist2        Distance to the centroid of that simplex
     */
    double dist2 = 0;  // squared distance
    for (int i = 0; i < D; i++) {
        dist2 += pow(target[i] - root->centroid[i], 2);
    }

    // Recursion exit condition - root is further than current closest neighbour
    if (dist2 - root->sbtree_radius_sq >= min_dist2) {
        return best;
    }

    // assert(root->lchild != NULL || root->rchild != NULL);

    // Perhaps the current point is the new closest?
    if (dist2 < min_dist2) {
        min_dist2 = dist2;
        best = root;
    }

    // Find closest child
    double ldist2 = 0;
    double rdist2 = 0;

    if (root->lchild != NULL) {
        for (int i = 0; i < D; i++) {
            ldist2 += pow(target[i] - root->lchild->centroid[i], 2);
        }
    }
    if (root->rchild != NULL) {
        for (int i = 0; i < D; i++) {
            rdist2 += pow(target[i] - root->rchild->centroid[i], 2);
        }
    }

    // Recurse into closest child first
    // I wonder if this would be prettier with gotos...
    Simplex<D> * old_best = best;
    if (ldist2 <= rdist2) {
        if (root->lchild != NULL) {
            best = find_nearest_neighbour_sbtree(root->lchild, target, best, min_dist2);
        }
        if (best != old_best) {
            // Recalculate min_dist2
            // TODO: Inefficient.
            min_dist2 = 0;
            for (int i = 0; i < D; i++) {
                min_dist2 += pow(target[i] - best->centroid[i], 2);
            }
//            std::cout << "New min_dist2: " << min_dist2 << std::endl;
        }
        if (root->rchild != NULL) {
            best = find_nearest_neighbour_sbtree(root->rchild, target, best, min_dist2);
        }
    } else {
        if (root->rchild != NULL) {
            best = find_nearest_neighbour_sbtree(root->rchild, target, best, min_dist2);
        }
        if (best != old_best) {
            // Recalculate min_dist2
            // TODO: Inefficient.
            min_dist2 = 0;
            for (int i = 0; i < D; i++) {
                min_dist2 += pow(target[i] - best->centroid[i], 2);
            }
        }
        if (root->lchild != NULL) {
            best = find_nearest_neighbour_sbtree(root->lchild, target, best, min_dist2);
        }
    }

    return best;
}


template<int N, int D, int S>
double Cool<N, D, S>::interpolate_sbtree(double * coords) {
    /**
     * Interpolate the given point.
     * First, find the closest simplex using the ball tree. Then, find the simplex containing the given point using
     * repeated "flips". Finally, interpolate using a weighted average (Delaunay).
     */

//    std::cout << "coords: ";
//    for (int i = 0; i < D; i++) {
//        std::cout << coords[i] << " ";
//    }
//    std::cout << std::endl;
//    std::map<Simplex<D> *, int> visited;
    Simplex<D> * best = NULL;
    Simplex<D> * nn = find_nearest_neighbour_sbtree(sbtree, coords, best, DBL_MAX);

    double * bary = nn->convert_to_bary(coords);
    int inside = nn->check_bary(bary);
    int flips = 0;

    while (!inside) {
        std::cout << flips << std::endl;
        // Calculate difference vectors from coords to midpoints; normalize; figure out the one with largest scalar
        // product (= smallest angle)
        double best_dir_dot = - DBL_MAX;
        int best_dir;
        for (int i = 0; i < D+1; i++) {
            double diff_vec[D];
            double diff_len = 0;
            for (int j = 0; j < D; j++) {
                diff_vec[j] = coords[j] - nn->midpoints[i][j];
                diff_len += pow(diff_vec[j], 2);
            }
            diff_len = sqrt(diff_len);
            for (int j = 0; j < D; j++) {
                diff_vec[j] /= diff_len;
            }

            double dot = 0;
            for (int j = 0; j < D; j++) {
                dot += diff_vec[j] * nn->normals[i][j];
            }

            if (dot > best_dir_dot) {
                best_dir_dot = dot;
                best_dir = i;
            }
        }

        nn = nn->neighbour_pointers[best_dir];
        std::cout << "Simplex index: ";
        int si = -1;
        for (int i = 0; i < S; i++) {
            if (&simplices[i] == nn) {
                si = i;
                break;
            }
        }
        std::cout << si << std::endl;
        bary = nn->convert_to_bary(coords);
        inside = nn->check_bary(bary);

        flips++;
        if (flips > 100) {
            std::cerr << "Error: More than 100 flips." << std::endl;
            abort();
            break;
        }

        std::cout << "inside: " << inside << std::endl;
//        // Error handling
//        if (best_dir_dot <= 0 && inside == 0) {
//            std::cerr << "Error in Simplex traversal: no positive scalar product: "
//                      << best_dir_dot << " <= 0" << std::endl;
//
//            std::cerr << "Simplex index: " << std::endl;
//            int si = -1;
//            for (int i = 0; i < S; i++) {
//                if (&simplices[i] == nn) {
//                    si = i;
//                    break;
//                }
//            }
//            if (si > -1) {
//                std::cerr << si << std::endl;
//                simplices[si].validate_normals();
//            } else {
//                std::cerr << "Not found" << std::endl;
//            }
//            std::cerr << "Target point in regular/barycentric coordinates:" << std::endl;
//            for (int i = 0; i < D; i++) {
//                std::cerr << coords[i] << " ";
//            }
//            std::cerr << "\t\t";
//            for (int i = 0; i < D; i++) {
//                std::cerr << bary[i] << " ";
//            }
//            std::cerr << std::endl;
//
//            std::cerr << "Simplex points:" << std::endl;
//            for (int j = 0; j < D+1; j++) {
//                std::cerr << j << ":\t";
//                for (int k = 0; k < D; k++) {
//                    std::cerr << nn->points[j]->coords[k] << " ";
//                }
//                std::cerr << std::endl;
//            }
//            std::cerr << "Centroid: " << std::endl;
//            for (int j = 0; j < D; j++) {
//                std::cerr << nn->centroid[j] << " ";
//            }
//            std::cerr << std::endl;
//            std::cerr << "Midpoints:" << std::endl;
//            for (int j = 0; j < D+1; j++) {
//                std::cerr << j << ":   ";
//                for (int k = 0; k < D; k++) {
//                    std::cerr << nn->midpoints[j][k] << " ";
//                }
//                std::cerr << std::endl;
//            }
//            std::cerr << "Midpoint diffs: " << std::endl;
//            for (int j = 0; j < D+1; j++) {
//                double diff[D];
//                double len = 0;
//                for (int k = 0; k < D; k++) {
//                    diff[k] = coords[k] - nn->midpoints[j][k];
//                    len += pow(diff[k], 2);
//                }
//                len = sqrt(len);
//                std::cerr << j << ":   ";
//                for (int k = 0; k < D; k++) {
//                    std::cerr << diff[k] / len << " ";
//                }
//                std::cerr << std::endl;
//            }
//            std::cerr << "Normal vectors:" << std::endl;
//            for (int j = 0; j < D+1; j++) {
//                std::cerr << j << ":\t";
//                for (int k = 0; k < D; k++) {
//                    std::cerr << nn->normals[j][k] << " ";
//                }
//                std::cerr << std::endl;
//            }
//
//            abort();
//        }
    }

//    int dbg_count = 0;
//    while (!inside) {
//        std::cout << dbg_count << std::endl;
//        for (int i = 0; i < S; i++) {
//            if (&simplices[i] == nn) {
//                std::cout << "Simplex: " << i << std::endl;
//                break;
//            }
//        }

//        std::cout << "Simplex points: " << std::endl;
//        for (int i = 0; i < D+1; i++) {
//            for (int j = 0; j < D; j++) {
//                std::cout << nn->points[i]->coords[j] << " ";
//            }
//            std::cout << std::endl;
//        }
//        std::cout << std::endl;

//        // Find nearest midpoint
//        std::cout << "Finding midpoints" << std::endl;
//        double midpoints[D+1][D];
//        for (int i = 0; i < D+1; i++) { // D+1 midpoints, one per plane
//            for (int j = 0; j < D; j++) {   // D coordinates
//                midpoints[i][j] = 0;
//                for (int k = 0; k < D; k++) {   // D points defining plane
//                    if (k < i) {
//                        midpoints[i][j] += nn->points[k]->coords[j];
//                    } else {
//                        midpoints[i][j] += nn->points[k+1]->coords[j];
//                    }
//                }
//                midpoints[i][j] /= D;
//            }
//        }
//        std::cout << "All midpoints: " << std::endl;
//        for (int i = 0; i < D+1; i++) {
//            for (int j = 0; j < D; j++) {
//                std::cout << midpoints[i][j] << " ";
//            }
//            std::cout << std::endl;
//        }

//        std::cout << "Projecting midpoint normals" << std::endl;
//        double best_dir_dot = 0;
//        int best_dir;
//        for (int i = 0; i < D+1; i++) {
//            double diff_vec[D];
//            double diff_len = 0;
//            for (int j = 0; j < D; j++) {
//                diff_vec[j] = coords[j] - midpoints[i][j];
//                diff_len += pow(diff_vec[j], 2);
//            }
//            diff_len = sqrt(diff_len);
//            for (int j = 0; j < D; j++) {
//                diff_vec[j] /= diff_len;
//            }
//
//            double dot = 0;
//            for (int j = 0; j < D; j++) {
//                dot += diff_vec[j] * nn->normals[i][j];
//            }
//            std::cout << "Midpoint: " << midpoints[i][0] << " " << midpoints[i][1] << std::endl;
//            std::cout << "Normal: " << nn->normals[i][0] << " " << nn->normals[i][1] << std::endl;
//            std::cout << "diff vec: " << diff_vec[0] << " " << diff_vec[1] << std::endl;
//            std::cout << std::endl;
//
//            if (dot > best_dir_dot) {
//                best_dir_dot = dot;
//                best_dir = i;
//            }
//            std::cout << "dot: " << dot << " best dot: " << best_dir_dot << std::endl;
//        }

//        std::cout << "Going to new nearest neighbor" << std::endl;
//        nn = nn->neighbour_pointers[best_dir];
//
//        std::cout << "convert to bary" << std::endl;
//        bary = nn->convert_to_bary(coords);
//
//        std::cout << "checking bary" << std::endl << std::endl;
//        inside = nn->check_bary(bary);
//
//        dbg_count++;
//        if (dbg_count > 100) {
//            std::cerr << "Error: More than 100 flips." << std::endl;
//            exit(1);
//            break;
//        }
//    }

//    std::cout << "Found simplex" << std::endl;

    // The actual interpolation step
    double val = 0;
    for (int i = 0; i < D+1; i++) {
        val += bary[i] * nn->points[i]->value;
    }

    delete[] bary;
    delete best;

    sbtree_interpolate_calls++;
    sbtree_flips += flips;
    avg_sbtree_flips = sbtree_flips / (float) sbtree_interpolate_calls;

    return val;
}


template<int N, int D, int S>
int Cool<N, D, S>::read_files(std::string cool_file, std::string tri_file, std::string neighbour_file) {
    /**
     * Reads files generated by the python program.
     *
     * cool_file        Path to file containing points.csv of the grid
     * tri_file         Path to file containing triangulation
     * neighbour_file    Path to file containing triangulation neighbourhood relations
     */
    std::ifstream file;
    std::string line;
    std::string value;


    /* Read points */
    file.open(cool_file);
    if (!file.is_open()) {
        std::cerr << "Error reading " << cool_file << std::endl;
        return 1;
    }

    for (int i = 0; i < N; i++) {
        std::getline(file, line);
        std::stringstream linestream(line);
        for (int j = 0; j < D; j++) {     // D coordinates, 1 value
            std::getline(linestream, value, ',');
            points[i].coords[j] = std::stod(value);
        }
        std::getline(linestream, value, ',');
        points[i].value = std::stod(value);

        for (int j = 0; j < D; j++) {
            if (points[i].coords[j] < mins[j]) {
//                std::cout << "Updating mins[" << j << "] from " << mins[j] << " to " << points[i].coords[j] << std::endl;
                mins[j] = points[i].coords[j];
            }
            if (points[i].coords[j] > maxs[j]) {
//                std::cout << "Updating maxs[" << j << "] from " << maxs[j] << " to " << points[i].coords[j] << std::endl;
                maxs[j] = points[i].coords[j];
            }
        }
    }


    file.close();

    /* Read Simplices */
    file.open(tri_file);
    if (!file.is_open()) {
        std::cerr << "Error reading " << tri_file << std::endl;
        return 2;
    }


    for (int i = 0; i < S; i++) {
        std::getline(file, line);
        std::stringstream linestream(line);

        int buffer[D+1];
        for (int j = 0; j < D+1; j++) {     // D+1 points per simplex
            std::getline(linestream, value, ',');
            buffer[j] = std::stoi(value);

            simplices[i].points[j] = &points[std::stoi(value)];
        }

        for (int j = 0; j < D+1; j++) {
            for (int k = 0; k < D+1; k++) {
                if (j != k && buffer[j] == buffer[k]) {
                    std::cerr << "Warning: Degenerate triangle (index " << i << ", points " << buffer[j] << " and " << buffer[k] << ")" << std::endl;
                }
            }
        }

        // scipy.spatial.Delaunay unfortunately likes making simplices out of the D-1d hyperplanes making up the
        // faces of the D dimensional hypercube which is the parameter space
        // But D-1d planes cannot be D dimensional simplices
        // This leads to problems when validating the normal vectors (and makes the normal vectors meaningless)
        // So those are skipped
//        simplices[i].validate_simplex();
        int skip = 0;
        for (int j = 0; j < D; j++) {
            int samemin = 1;
            int samemax = 1;
            for (int k = 0; k < D+1; k++) {
                if (simplices[i].points[k]->coords[j] != mins[j]) {
                    samemin = 0;
                    break;
                }
            }
            for (int k = 0; k < D+1; k++) {
                if (simplices[i].points[k]->coords[j] != maxs[j]) {
                    samemax = 0;
                    break;
                }
            }
            if (samemin || samemax) {
                skip = 1;
                break;
            }
        }

        simplices[i].calculate_centroid();
        simplices[i].calculate_midpoints();

        for (int j = 0; j < D; j++) {
            std::cout << mins[j] << " " << maxs[j] << std::endl;
        }
        if (skip == 1) {
            continue;
        }
        simplices[i].calculate_normals();
        simplices[i].validate_normals();
    }


    file.close();

    // Construct matrices for each simplex
    // https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Barycentric_coordinates_on_tetrahedra
    for (int s = 0; s < S; s++) {   // for each simplex

        for (int i = 0; i < D; i++) {   // for each point (except the last)
            for (int j = 0; j < D; j++) {   // for each coordinate
                // matrix[j][i], not matrix[i][j] - coordinates go down, points go right
                simplices[s].T_inv[j][i] = simplices[s].points[i]->coords[j] - simplices[s].points[D]->coords[j];
            }
        }

        // Right now its just T - now invert
        simplices[s].invert_T();

    }

    /* Read neighbourhood relations */
    file.open(neighbour_file);
    if (!file.is_open()) {
        std::cerr << "Error reading " << neighbour_file << std::endl;
        return 2;
    }


    for (int i = 0; i < S; i++) {
        std::getline(file, line);
        std::stringstream linestream(line);
        for (int j = 0; j < D+1; j++) {     // D+1 points.csv per simplex
            std::getline(linestream, value, ',');
            simplices[i].neighbour_indices[j] = std::stoi(value);

            if (std::stoi(value) != -1) {
                simplices[i].neighbour_pointers[j] = &(simplices[std::stoi(value)]);
            }
        }
    }

    file.close();

//    // Shrink point vectors
//    for (int i = 0; i < N; i++) {
//        points[i].simplices.shrink_to_fit();
//    }
    return 0;
}


template<int D>
void Simplex<D>::invert_T() {
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

    for (int i = 0; i < D; i++) {
        for (int j = 0; j < D; j++) {
            T_inv[i][j] = unit[i][j];
        }
    }
}



template<int D>
template<size_t M, size_t N1, size_t N2>
void Simplex<D>::gauss_elimination(double (&A)[M][N1], double (&B)[M][N2], int apply_permutation) {
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
     *                          0 for now, 1 for yes
     */

    // Permutation of ROW InDiXes
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
        if (pivot == 0) {
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

        for (int j = 0; j < M; j++) {
            for (int k = 0; k < N1; k++) {
                std::cout << A[j][k] << " ";
            }
            std::cout << "| ";
            for (int k = 0; k < N2; k++) {
                std::cout << B[j][k] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    if (apply_permutation) {
        for (int i = 0; i < M; i++) {
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
        }
    }
}


//template<int D>
//double * Simplex<D>::gauss_elimination(double * A, double * B, int M, int N) {
//    /**
//     * Ax = B
//     * Overloaded version of gauss_elimination, for when B is only a vector.
//     */
//     return gauss_elimination(A, B, M, N, 1);
//}


template<int D>
double * Simplex<D>::convert_to_bary(const double * coords) {
    /**
     * Convert the given coordinates to the barycentric coordinate system of the simplex.
     *
     * https://math.stackexchange.com/a/1226825
     */
    // 1. Obtain p - v_n+1
    double vec[D];
    for (int i = 0; i < D; i++) {
        vec[i] = coords[i] - points[D]->coords[i];
    }

    // 2. Multiply T_inv * (p - v_n+1) to get lambda
    double * bary = new double[D+1];
    for (int i = 0; i < D; i++) {
        bary[i] = 0;
        for (int j = 0; j < D; j++) {
            bary[i] += T_inv[i][j] * vec[j];
        }
    }

    // dependent coordinate lambda_n+1
    bary[D] = 1;
    for (int i = 0; i < D; i++) {
        bary[D] -= bary[i];
    }

    return bary;
}


template<int D>
inline int Simplex<D>::check_bary(const double* bary) {
    /**
     * Check if the given barycentric coordinates belong to a point inside or outside of the simplex
     * Returns 1 if inside, 0 otherwise.
     * Edge cases are considered to be outside. TODO avoids edge cases but means points.csv are never associated with anything?
     */
    double bsum = 0;
    for (int i = 0; i < D+1; i++) {
        if (bary[i] < 0) {
            std::cout << i << ": bary " << bary[i] << " < 0" << std::endl;
            // Consider points on edges to not be contained, in order to avoid degeneracies!
            return 0;
        }
        if (isnanf(bary[i])) {
            std::cerr << "Error: Barycentric coordinate " << i << " is nan" << std::endl;
//            return 0;
            abort();
        }
        bsum += bary[i];
        std::cout << bsum << std::endl;
    }

    if (bsum > 1 + EPSILON) {
        for (int i = 0; i < D+1; i++) {
            std::cout << bary[i] << " + ";
        }
        std::cout << " = " << bsum << " > 1" << std::endl;
        return 0;
    }

    return 1;
}


template<int N>
class MultilinearInterpolator {
private:
    int dims[N];
    double * grid;
    double minmax[N][2];
    double dim_lens[N];
    double diffs[N];
    double dim_offsets[N];

    double linear_interpolation(double, double, double, double, double);
public:
    MultilinearInterpolator(double * in_grid, const double * in_minmax, const int * in_dims) {
        grid = in_grid;
        for (int i = 0; i < N; i++) {
            minmax[i][0] = in_minmax[2*i];
            minmax[i][1] = in_minmax[2*i+1];
            dim_lens[i] = minmax[i][1] - minmax[i][0];
        }
        for (int i = 0; i < N; i++) {
            dims[i] = in_dims[i];
            diffs[i] = dims[i] / dim_lens[i];
        }
        for (int i = 0; i < N; i++) {
            dim_offsets[i] = 1;
            if (i == N) {
                break;
            }
            for (int j = i+1; j < N; j++) {
                dim_offsets[i] *= dims[j];
            }
        }
    }
    double interpolate(const double *);
};


template<int N>
double MultilinearInterpolator<N>::linear_interpolation(double x, double x0, double y0, double x1, double y1) {
    return (y0 * (x1 - x) + y1 * (x - x0)) / (x1 - x0);
}

template<int N>
double MultilinearInterpolator<N>::interpolate(const double* point) {
    /**
     * Optimized multilinear interpolation function.
     * Undefined behaviour when trying to interpolate outside grid.
     *
     * point = array of length N, point to interpolate at
     */

    // 1. Figure out points of surrounding hypercube
    int indices[N];

    for (int i = 0; i < N; i++) {
        // Automatically cast to int!
        indices[i] = (point[i] - minmax[i][0]) * dims[i] / dim_lens[i];
    }

    double * vals = new double[(int) pow(2, N)];
    for (int i = 0; i < pow(2, N); i++) {
        std::bitset<N> offsets(i);
        int index = 0;

        for (int j = 0; j < N; j++) {
            // TODO: Fix
            index += (offsets[j]) ? dim_offsets[j] * (indices[j] + 1) : dim_offsets[j] * indices[j];
        }

        vals[i] = grid[index];
        std::cout << offsets << " " << vals[i] << std::endl;
    }

    double * new_vals;
    for (int i = N; i > 0; i--) {
        new_vals = new double[(int) pow(2, i-1)];
        for (int j = 0; j < (int) pow(2, i-1); j++) {
            // To understand why the offset in the second vals[] works, look at the bit representation of the indices
            // Those bits correspond to the offsets above
            new_vals[j] = vals[j] * (1 - diffs[i]) + vals[j + (int) pow(2, i-1)] * diffs[i];
        }

        delete[] vals;
        vals = new_vals;
    }

    double result = new_vals[0];
    delete[] new_vals;
    return result;
}
