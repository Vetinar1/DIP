//
// Created by vetinari on 26.08.20.
//

#include "main.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <float.h>

int main() {
    double * points     = (double*)malloc(n_points * (cool_dim+1) * sizeof(double));
//    int * triangulation = (int*)   malloc(n_simplices * (cool_dim+1) * sizeof(int));
    int st_s = sizeof(simplex_t);
    // Pointer to "array" of pointers to simplices
    simplex_t ** triangulation = (simplex_t**)malloc(n_simplices * sizeof(simplex_t *));

    // TODO: File reading as generic function (dynamic input size)
    // TODO: Do not require file line number as input. realloc?
    FILE * f = fopen("../data.csv", "r");

    // Read points
    if (f == NULL) {
        perror("Error opening file data.csv: ");
        return 1;
    }
    for (int i = 0; i < n_points; i++) {
        fscanf(f, "%lf,%lf,%lf\n", points + i*(cool_dim+1), points + i*(cool_dim+1) + 1, points + i*(cool_dim+1) + 2);
    }

    for (int i = 0; i < n_points; i++) {
        printf("%f %f %f\n", *(points + i*(cool_dim+1)), *(points + i*(cool_dim+1) + 1), *(points + i*(cool_dim+1) + 2));
    }

    int close = fclose(f);
    if (close != 0) {
        perror("Error closing file data.csv: ");
        return 1;
    }

    // Read Triangulation
    f = fopen("../dtriangulation", "r");
    if (f == NULL) {
        perror("Error opening file dtriangulation: ");
        return 1;
    }
    for (int i = 0; i < n_points; i++) {
        simplex_t * new_simplex = (simplex_t*)malloc(sizeof(simplex_t)); // TODO free
        fscanf(f, "%d %d %d\n", new_simplex->points, new_simplex->points + 1, new_simplex->points + 2);

        for (int j = 0; j < cool_dim; j++) {
            new_simplex->centroid[j] = 0;
            new_simplex->part_of_tree = 0;

            for (int k = 0; k < cool_dim+1; k++) {
                new_simplex->centroid[j] += *(points + (cool_dim+1) * new_simplex->points[k] + j);
            }

            new_simplex->centroid[j] /= (cool_dim+1);
            *(triangulation + i*sizeof(simplex_t*)) = new_simplex;
        }
    }

    simplex_t * root = build_ball_tree(n_simplices, triangulation);

    return 0;

    int N = 2;
    int count = 1;
    double matrix[2][2] = {
            {-1.25, 1.25},
            {-2., -1.5}
    };

    invert_matrix(&matrix[0][0], N);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
    return 0;
}

simplex_t * build_ball_tree(int N, simplex_t ** triangulation) {
    /*
     * https://en.wikipedia.org/wiki/Ball_tree#k-d_Construction_Algorithm
     * N                number of simplices left
     * triangulation    Array of pointers to simplices
     *
     * A lot of the issues in this function are redundant and could be merged, but I don't think it will be an issue.
     */

    // Single point?
    assert(N > 0);
    if (N == 1) {
        (*(triangulation))->btree_radius_sq = 0;
        (*(triangulation))->lchild = NULL;
        (*(triangulation))->rchild = NULL;
        return (*(triangulation));
    }

    // N > 1
    // Find dimension of greatest spread
    double largest_spread = -1;
    int spread_dim = -1;
    double spread_dim_min = -1;
    for (int i = 0; i < cool_dim; i++) {
        double dim_min = 0, dim_max = 0, spread = 0;
        for (int j = 0; j < N; j++) {
            double val = (*(triangulation + j * sizeof(simplex_t*)))->centroid[i];
            if (val > dim_max) {
                dim_max = val;
            } else if (val < dim_min) {
                dim_min = val;
            }
        }

        spread = dim_max - dim_min;
        if (spread > largest_spread) {
            largest_spread = spread;
            spread_dim = i;
            spread_dim_min = dim_min;
        }
    }
    assert(largest_spread != -1 && spread_dim != -1 && spread_dim_min != -1);

    // Find midpoint - point closest to average, i.e. dim_min + 0.5 * spread
    double avg = spread_dim_min + 0.5 * largest_spread;
    double min_dist = DBL_MAX;
    int min_dist_index = -1;
    for (int i = 0; i < N; i++) {
        double val = (*(triangulation + i * sizeof(simplex_t*)))->centroid[spread_dim];
        double dist = fabs(avg - val);

        if (dist < min_dist) {
            min_dist = dist;
            min_dist_index = i;
        }
    }
    assert(min_dist_index != -1 && min_dist != DBL_MAX);
    simplex_t * pivot_addr = *(triangulation + min_dist_index * sizeof(simplex_t*));

    // Find sets of simplices to the left and right of the midpoint along spread dimension
    simplex_t ** L = (simplex_t **)malloc(N*sizeof(simplex_t*));
    simplex_t ** R = (simplex_t **)malloc(N*sizeof(simplex_t*));
    int lcount = 0;
    int rcount = 0;

    for (int i = 0; i < N; i++) {
        simplex_t * addr = *(triangulation + i*sizeof(simplex_t*));
        if (addr->centroid[spread_dim] < avg) {
            *(L + lcount * sizeof(simplex_t*)) = addr;
            lcount++;
        } else if (addr->centroid[spread_dim] > avg) {
            *(R + rcount * sizeof(simplex_t*)) = addr;
            rcount++;
        }
    }

    // Recursion
    if (lcount > 0) {
        pivot_addr->lchild = build_ball_tree(lcount, L);
    } else {
        pivot_addr->lchild = NULL;
    }
    if (rcount > 0) {
        pivot_addr->rchild = build_ball_tree(rcount, R);
    } else {
        pivot_addr->rchild = NULL;
    }

    // Find (square of) radius of the ball of current element
    // Note that unlike min_dist, which is only in the spread dimension, this is in the full space!
    double max_dist = 0;
    for (int i = 0; i < N; i++) {
        double dist = 0;

        for (int j = 0; j < cool_dim; j++) {
//            dist += pow(pivot_addr->centroid[j] - (*(triangulation + i*sizeof(simplex_t*)))->centroid[i], 2);
            dist += (pivot_addr->centroid[j] - (*(triangulation + i*sizeof(simplex_t*)))->centroid[i]) * (pivot_addr->centroid[j] - (*(triangulation + i*sizeof(simplex_t*)))->centroid[i]);
        }

        if (dist > max_dist) {
            max_dist = dist;
        }
    }

    assert(max_dist > 0);

    pivot_addr->btree_radius_sq = max_dist;

    return pivot_addr;
}

void invert_matrix(double* matrix, int N) {
    /*
     * Invert Matrix in place using Gauss Jordan Elimination with partial pivoting
     *
     * Algorithm adapted from Numerical Recipes, 2.1
     *
     * TODO frequently called, better be efficient
     * TODO implicit pivoting since we are not in a normalized space
     */

    int row_idx [N];
    for (int i = 0; i < N; i++) {
        row_idx[i] = i;
    }

    // yes, you can work in place instead, but i cant be bothered right now
    double * unit = calloc(N * N, sizeof(double));
    for (int i = 0; i < N; i++) {
        *(unit + N * i + i) = 1;
    }


    for (int i = 0; i < N; i++) {
        // find pivot
        int pivot_row_idx = i;
        double pivot = matrix[N*row_idx[i] + i];
        for (int j = i; j < N; j++) {
            if (fabs(matrix[N*row_idx[j]+i]) > fabs(pivot)) {
                pivot = matrix[N*row_idx[j]+i];
                pivot_row_idx = j;
            }
        }

        // pivoting
        int temp = row_idx[i];
        row_idx[i] = row_idx[pivot_row_idx];
        row_idx[pivot_row_idx] = temp;

        // Normalization of current row
        double norm_factor = matrix[N*row_idx[i] + i];
        for (int j = 0; j < N; j++) {
            matrix[N*row_idx[i] + j] /= norm_factor;
            unit[N*row_idx[i] + j]   /= norm_factor;
        }

        // Subtracting the current row from all other rows
        for (int j = 0; j < N; j++) {
            if (j != i) {
                double factor = matrix[N * row_idx[j] + i];
                for (int k = 0; k < N; k++) {
                    matrix[N * row_idx[j] + k] -= factor * matrix[N * row_idx[i] + k];
                    unit[N * row_idx[j] + k] -= factor * unit[N * row_idx[i] + k];
                }
            }
        }

//        for (int j = 0; j < N; j++) {
//            for (int k = 0; k < N; k++) {
//                printf("%f ", matrix[N*row_idx[j] + k]);
//            }
//            printf("\n");
//        }
//        printf("\n");
//        for (int j = 0; j < N; j++) {
//            for (int k = 0; k < N; k++) {
//                printf("%f ", unit[N*row_idx[j] + k]);
//            }
//            printf("\n");
//        }
//        printf("\n\n\n");
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
//            matrix[N*row_idx[i]+j] = unit[N*row_idx[i]+j];
            matrix[N*row_idx[i]+j] = unit[N*i+j]; // TODO verify
//            matrix[N*i+j] = unit[N*i+j];
        }
    }
}

void solve_linear(double* matrix, double* vector, int N) {
    /*
     * Invert Matrix in place using Gauss Jordan Elimination with partial pivoting
     *
     * Algorithm adapted from Numerical Recipes, 2.1
     *
     * TODO frequently called, better be efficient
     * TODO implicit pivoting since we are not in a normalized space
     *
     * TODO TEST
     */

    int row_idx [N];
    for (int i = 0; i < N; i++) {
        row_idx[i] = i;
    }

    for (int i = 0; i < N; i++) {
        // find pivot
        int pivot_row_idx = i;
        double pivot = matrix[N*row_idx[i] + i];
        for (int j = i; j < N; j++) {
            if (fabs(matrix[N*row_idx[j]+i]) > fabs(pivot)) {
                pivot = matrix[N*row_idx[j]+i];
                pivot_row_idx = j;
            }
        }

        // pivoting
        int temp = row_idx[i];
        row_idx[i] = row_idx[pivot_row_idx];
        row_idx[pivot_row_idx] = temp;

        // Normalization of current row
        double norm_factor = matrix[N*row_idx[i] + i];
        for (int j = 0; j < N; j++) {
            matrix[N*row_idx[i] + j] /= norm_factor;
            vector[N*row_idx[i] + j] /= norm_factor;
        }

        // Subtracting the current row from all other rows
        for (int j = 0; j < N; j++) {
            if (j != i) {
                double factor = matrix[N * row_idx[j] + i];
                for (int k = 0; k < N; k++) {
                    matrix[N * row_idx[j] + k] -= factor * matrix[N * row_idx[i] + k];
                    vector[N * row_idx[j] + k] -= factor * vector[N * row_idx[i] + k];
                }
            }
        }

//        for (int j = 0; j < N; j++) {
//            for (int k = 0; k < N; k++) {
//                printf("%f ", matrix[N*row_idx[j] + k]);
//            }
//            printf("\n");
//        }
//        printf("\n");
//        for (int j = 0; j < N; j++) {
//            for (int k = 0; k < N; k++) {
//                printf("%f ", unit[N*row_idx[j] + k]);
//            }
//            printf("\n");
//        }
//        printf("\n\n\n");
    }
}


int contains(simplex_t * tri, const double * points, int pt_idx) {
    /*
     * Check if point in points at index pt_index is contained by simplex_t simplex.
     * Return 1 if true.
     * Return -1 if false.
     *
     * https://math.stackexchange.com/a/1226825
     */

    // 1. Build Matrix T to invert
    double matrix[cool_dim][cool_dim];
    for (int i = 0; i < cool_dim; i++) {
        // For each row
        for (int j = 0; j < cool_dim; j++) {
            // For each column
            matrix[i][j] = *(points + (cool_dim+1)*tri->points[j] + i) -
                           *(points + (cool_dim+1)*tri->points[cool_dim] + i);
        }
    }
    printf("\nMatrix T:\n");
    for (int i = 0; i < cool_dim; i++) {
        for (int j = 0; j < cool_dim; j++) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }

    // 2. Invert T
    invert_matrix(&matrix[0][0], cool_dim);
    printf("\nInverse T:\n");
    for (int i = 0; i < cool_dim; i++) {
        for (int j = 0; j < cool_dim; j++) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }

    // 3. Obtain p - v_n+1
    double vec[cool_dim];
    for (int i = 0; i < cool_dim; i++) {
        vec[i] = *(points + (cool_dim+1)*pt_idx + i) - *(points + (cool_dim+1)*tri->points[cool_dim] + i);
    }

    // 4. Multiply T_inv * (p - v_n+1) to get lambda
    double bary[cool_dim+1];
    for (int i = 0; i < cool_dim; i++) {
        bary[i] = 0;
        for (int j = 0; j < cool_dim; j++) {
            bary[i] += matrix[i][j] * vec[j];
        }
    }

    // dependent coordinate lambda_n+1
    bary[cool_dim] = 1;
    for (int i = 0; i < cool_dim; i++) {
        bary[cool_dim] -= bary[i];
    }

    // Check barycentric coordinates to see if point is contained
    printf("Barycentric coordinates: ");
    for (int i = 0; i < cool_dim+1; i++) {
        printf("%f ", bary[i]);
    }
    printf("\n");

    double bsum = 0;
    for (int i = 0; i < cool_dim+1; i++) {
        if (bary[i] <= 0) {
            // Consider points on edges to not be contained, in order to avoid degeneracies!
            printf("Not contained; Coordinate below zero\n");
            return -1;
        }
        bsum += bary[i];
    }

    if (bsum > 1) {
        printf("Not contained; Sum of coordinates greater than 1\n");
        return -1;
    }

    printf("Point contained\n");
    return 1;
}


//
//// Note: I'm keeping these outdated functions around in case they become useful at some point
//double signed_dist_to_facet(qfacet_t * qf, double * p, int d) {
//    /*
//     * calculate signed distance from facet struct to point found at pointer p. in space of dimension d
//     */
//    // Note that plane offset needs to be taken into account. (Draw a diagram)
//    // Since its a plane, ANY point will work as offset, so i just use the first of the facet
//    // TODO: Verify
//    double * subtr;
//    for (int k = 0; k < d; k++) {
//        subtr[k] = p[k] - qf->points[0];
//    }
//
//    double signed_dist = 0;
//    for (int k = 0; k < d; k++) {
//        signed_dist += subtr[k] * qf->normal[k];
//    }
//
//    return signed_dist;
//}
//
//
//// TODO N and d mean different things here than in the other function, unify
//void calculate_normal(qfacet_t * qf, double * points, int N, int d) {
//    /*
//     * Calculate normal of a qfacet
//     * https://math.stackexchange.com/a/3398303
//     *
//     */
//    // Build matrix
//    double * matrix = (double*)calloc(N*N, sizeof(double));
//
//    for (int i = 0; i < N; i++) {
//        for (int j = 0; j < N; j++) {
//            *(matrix + N*i + j) = *(points + d*qf->points[i] + j);
//        }
//
//        qf->normal[i] = 1;
//    }
//
//    solve_linear(&matrix[0], &qf->normal[0], N);
//
//    double length = 0;
//    for (int i = 0; i < N; i++) {
//        length += qf->normal[i] * qf->normal[i];
//    }
//    length = sqrt(length);
//    for (int i = 0; i < N; i++) {
//        qf->normal[i] /= length;
//    }
//}


