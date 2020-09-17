//
// Created by vetinari on 26.08.20.
//

#include "main.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <string.h>


int main() {
    double * points     = (double*)malloc(n_points * (cool_dim+1) * sizeof(double));
    int st_s = sizeof(simplex_t);
    // Pointer to "array" of pointers to simplices
    simplex_t ** triangulation = (simplex_t**)malloc(n_simplices * sizeof(simplex_t *));

    // TODO: File reading as generic function (dynamic input size)
    // TODO: Do not require file line number as input. realloc?
    FILE * f = fopen("../points", "r");

    // Read points
    if (f == NULL) {
        perror("Error opening file points: ");
        return 1;
    }
    for (int i = 0; i < n_points; i++) {
//        fscanf(f, "%lf %lf %lf", points + (i*(cool_dim+1))*sizeof(double), points + (i*(cool_dim+1) + 1)*sizeof(double), points + (i*(cool_dim+1) + 2)*sizeof(double));
        double temp1, temp2, temp3;
        fscanf(f, "%lf %lf %lf", &temp1, &temp2, &temp3);
        printf("Read: %i %f %f %f\n", i, temp1, temp2, temp3);
        *(points + (i*(cool_dim+1))*sizeof(double)) = temp1;
        *(points + (i*(cool_dim+1) + 1)*sizeof(double)) = temp2;
        *(points + (i*(cool_dim+1) + 2)*sizeof(double)) = temp3;
    }

    for (int i = 0; i < n_points; i++) {
        printf("%i %lf %lf %lf \n", i, *(points + ((cool_dim+1)*i)*sizeof(double)), *(points + ((cool_dim+1)*i+1)*sizeof(double)), *(points + ((cool_dim+1)*i+2)*sizeof(double)));
    }

    int close = fclose(f);
    if (close != 0) {
        perror("Error closing file points: ");
        return 1;
    }


    printf("2\n");
    // Read Triangulation
    f = fopen("../dtriangulation", "r");
    if (f == NULL) {
        perror("Error opening file dtriangulation: ");
        return 1;
    }

    for (int i = 0; i < n_simplices; i++) {
        simplex_t * new_simplex = (simplex_t*)malloc(sizeof(simplex_t));
        fscanf(f, "%d %d %d\n", new_simplex->points, new_simplex->points + 1, new_simplex->points + 2);

        for (int j = 0; j < cool_dim; j++) {
            new_simplex->centroid[j] = 0;
            for (int k = 0; k < cool_dim+1; k++) {
                new_simplex->centroid[j] += *(points + ((cool_dim+1) * new_simplex->points[k] + j)*sizeof(double));
            }

            new_simplex->centroid[j] /= (cool_dim+1);
            *(triangulation + i*sizeof(simplex_t*)) = new_simplex;
        }

    }


    simplex_t ** tricopy = malloc(n_simplices*sizeof(simplex_t*));
    memcpy(tricopy, triangulation, n_simplices*sizeof(simplex_t*));
    for (int i = 0; i < n_simplices; i++) {
        printf("%i\n", i);
        assert(*(tricopy + i*sizeof(simplex_t*)) == *(triangulation + i*sizeof(simplex_t*)));
    }
    return 0;

    simplex_t * root = build_ball_tree_debug(n_simplices, triangulation, 0);

    return 0;


    // Old, ignore:
//    int N = 2;
//    int count = 1;
//    double matrix[2][2] = {
//            {-1.25, 1.25},
//            {-2., -1.5}
//    };
//
//    invert_matrix(&matrix[0][0], N);
//
//    for (int i = 0; i < N; i++) {
//        for (int j = 0; j < N; j++) {
//            printf("%f ", matrix[i][j]);
//        }
//        printf("\n");
//    }
//    return 0;
}

simplex_t * build_ball_tree_debug(int N, simplex_t ** triangulation, int depth) {
    depth++;

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
    int spread_dim = 1;
    double avg = 0.1;
    double min_dist = DBL_MAX;
    int pivot_index = -1;
    for (int i = 0; i < N; i++) {
        double val = (*(triangulation + i * sizeof(simplex_t*)))->centroid[spread_dim];
        double dist = fabs(avg - val);

        if (dist < min_dist) {
            min_dist = dist;
            pivot_index = i;
//            printf("New midpoint: %i\n", pivot_index);
        }
    }
    assert(pivot_index != -1 && min_dist != DBL_MAX);
    simplex_t * pivot_addr = *(triangulation + pivot_index * sizeof(simplex_t*));

    // Find sets of simplices to the left and right of the midpoint along spread dimension
    simplex_t ** tcopy = (simplex_t**)malloc(N*sizeof(simplex_t*));
    memcpy(tcopy, triangulation, N*sizeof(simplex_t*));
    simplex_t ** L = (simplex_t **)malloc(N*sizeof(simplex_t*));
    simplex_t ** R = (simplex_t **)malloc(N*sizeof(simplex_t*));
    int lcount = 0;
    int rcount = 0;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i != j) {
//                printf("depth %i i %i j %i\n", depth, i, j);
                assert(*(triangulation + i*sizeof(simplex_t*)) != *(triangulation + j*sizeof(simplex_t*)));
            }
        }
    }
    for (int i = 0; i < N; i++) {
        printf("%i %i\n", depth, i);
        assert(*(triangulation + i*sizeof(simplex_t*)) == *(tcopy + i*sizeof(simplex_t*)));
    }

    for (int i = 0; i < N; i++) {
        simplex_t * addr = *(triangulation + i*sizeof(simplex_t*));

//        for (int j = 0; j <N; j++) {
//            printf("%i %p ", j, (void*)*(triangulation + j*sizeof(simplex_t*)));
//        }
//        printf("\n");
//        printf("%i %i %i %i %p\n", depth, i, lcount, rcount, (void*)*(triangulation + 106*sizeof(simplex_t*)));
        for (int j = 0; j < N; j++) {
            if (i != j) {
//                printf("%i %i %i %p %p\n", depth, i, j, (void*)addr, (void*)*(triangulation + j*sizeof(simplex_t*)));
//                printf("%i %i %i \n", depth, i, j );
                if (addr == *(triangulation + j*sizeof(simplex_t*))) {
                    printf("222 %i %i\n", i, j);
                }
                assert(addr != *(triangulation + j*sizeof(simplex_t*)));
            } else {
                assert(addr == *(triangulation + i*sizeof(simplex_t*)));
            }
        }

        if (i == pivot_index)
            continue;

        if (addr->centroid[spread_dim] <= avg) {    // If in doubt left side - TODO: Problem?
            *(L + lcount * sizeof(simplex_t*)) = addr;
            lcount++;
        } else if (addr->centroid[spread_dim] > avg) {
            *(R + rcount * sizeof(simplex_t*)) = addr;
            rcount++;
        }
    }

    printf("lc %i rc %i N %i\n", lcount, rcount, N);
    assert(lcount + rcount + 1 == N);
//    printf("lcount: %i , rcount: %i\n", lcount, rcount);

    // Recursion
    if (lcount > 0) {
        printf("Recursing at %i into L (%i)\n", depth, lcount);
        pivot_addr->lchild = build_ball_tree_debug(lcount, L, depth);
    } else {
        pivot_addr->lchild = NULL;
    }
    if (rcount > 0) {
        printf("Recursing at %i into R (%i)\n", depth, rcount);
        pivot_addr->rchild = build_ball_tree_debug(rcount, R, depth);
    } else {
        pivot_addr->rchild = NULL;
    }
    printf("Finished at %i\n", depth);

    return pivot_addr;
}

simplex_t * build_ball_tree(int N, simplex_t ** triangulation, int depth) {
    /*
     * https://en.wikipedia.org/wiki/Ball_tree#k-d_Construction_Algorithm
     * N                number of simplices left
     * triangulation    Array of pointers to simplices
     *
     * A lot of the issues in this function are redundant and could be merged, but I don't think it will be an issue.
     */
    depth++;

//    for (int i = 0; i < N; i++) {
//        printf("%p\n", *(triangulation + i*sizeof(simplex_t*)));
//    }

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
    double spread_dim_min = DBL_MAX;
    double spread_dim_max = -1*DBL_MAX;
    for (int i = 0; i < cool_dim; i++) {
        double dim_min = DBL_MAX, dim_max = -1*DBL_MAX, spread = 0;
        for (int j = 0; j < N; j++) {
            double val = (*(triangulation + j * sizeof(simplex_t*)))->centroid[i];
//            printf("val %i : %lf\n", j, val);
            if (val > dim_max) {
                dim_max = val;
            }
            if (val < dim_min) {
                dim_min = val;
            }
        }

        spread = dim_max - dim_min;
        if (spread > largest_spread) {
            largest_spread = spread;
            spread_dim = i;
            spread_dim_min = dim_min;
            spread_dim_max = dim_max;
        }
    }
    assert(largest_spread != -1 && spread_dim != -1 && spread_dim_min != -1);

//    printf("Largest spread: %lf, dimension: %i\n", largest_spread, spread_dim);
//    printf("Largest spread: %lf\n", largest_spread);
//    printf("Spread dim: %i\n", spread_dim);
//    printf("Spread dim min: %lf\n", spread_dim_min);
//    printf("Spread dim max: %lf\n", spread_dim_max);

    // Find midpoint - point closest to average, i.e. dim_min + 0.5 * spread
    double avg = spread_dim_min + 0.5 * largest_spread;
//    printf("avg %lf \n", avg);
    double min_dist = DBL_MAX;
    int pivot_index = -1;
    for (int i = 0; i < N; i++) {
        double val = (*(triangulation + i * sizeof(simplex_t*)))->centroid[spread_dim];
        double dist = fabs(avg - val);

        if (dist < min_dist) {
            min_dist = dist;
            pivot_index = i;
//            printf("New midpoint: %i\n", pivot_index);
        }
    }
    assert(pivot_index != -1 && min_dist != DBL_MAX);
    simplex_t * pivot_addr = *(triangulation + pivot_index * sizeof(simplex_t*));

    // Find sets of simplices to the left and right of the midpoint along spread dimension
    simplex_t ** L = (simplex_t **)malloc(N*sizeof(simplex_t*));
    simplex_t ** R = (simplex_t **)malloc(N*sizeof(simplex_t*));
    int lcount = 0;
    int rcount = 0;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i != j) {
//                printf("depth %i i %i j %i\n", depth, i, j);
                assert(*(triangulation + i*sizeof(simplex_t*)) != *(triangulation + j*sizeof(simplex_t*)));
            }
        }
    }

    for (int i = 0; i < N; i++) {
        simplex_t * addr = *(triangulation + i*sizeof(simplex_t*));

//        for (int j = 0; j <N; j++) {
//            printf("%i %p ", j, (void*)*(triangulation + j*sizeof(simplex_t*)));
//        }
//        printf("\n");
        printf("%i %i %i %i %p\n", depth, i, lcount, rcount, (void*)*(triangulation + 106*sizeof(simplex_t*)));
        for (int j = 0; j < N; j++) {
            if (i != j) {
//                printf("%i %i %i %p %p\n", depth, i, j, (void*)addr, (void*)*(triangulation + j*sizeof(simplex_t*)));
//                printf("%i %i %i \n", depth, i, j );
                if (addr == *(triangulation + j*sizeof(simplex_t*))) {
                    printf("%i %i\n", i, j);
                }
                assert(addr != *(triangulation + j*sizeof(simplex_t*)));
            } else {
                assert(addr == *(triangulation + i*sizeof(simplex_t*)));
            }
        }
        if (i == pivot_index)
            continue;

        if (addr->centroid[spread_dim] <= avg) {    // If in doubt left side - TODO: Problem?
            *(L + lcount * sizeof(simplex_t*)) = addr;
            lcount++;
        } else if (addr->centroid[spread_dim] > avg) {
            *(R + rcount * sizeof(simplex_t*)) = addr;
            rcount++;
        }
    }

    for (int i = 0; i < lcount; i++) {
        for (int j = 0; j < lcount; j++) {
            if (i != j) {
//                printf("depth %i i %i j %i\n", depth, i, j);
                assert(*(L + i*sizeof(simplex_t*)) != *(L + j*sizeof(simplex_t*)));
            }
        }
    }
    for (int i = 0; i < rcount; i++) {
        for (int j = 0; j < rcount; j++) {
            if (i != j) {
//                printf("depth %i i %i j %i\n", depth, i, j);
                assert(*(R + i*sizeof(simplex_t*)) != *(R + j*sizeof(simplex_t*)));
            }
        }
    }
    printf("lc %i rc %i N %i\n", lcount, rcount, N);
    assert(lcount + rcount + 1 == N);
//    printf("lcount: %i , rcount: %i\n", lcount, rcount);

    // Recursion
    if (lcount > 0) {
        printf("Recursing at %i into L (%i)\n", depth, lcount);
        pivot_addr->lchild = build_ball_tree(lcount, L, depth);
    } else {
        pivot_addr->lchild = NULL;
    }
    if (rcount > 0) {
        printf("Recursing at %i into R (%i)\n", depth, rcount);
        pivot_addr->rchild = build_ball_tree(rcount, R, depth);
    } else {
        pivot_addr->rchild = NULL;
    }
    printf("Finished at %i\n", depth);

    // Find (square of) radius of the ball of current element
    // Note that unlike min_dist, which is only in the spread dimension, this is in the full space!
    double max_dist = 0;
    for (int i = 0; i < N; i++) {
        if (i == pivot_index)
            continue;

        double dist = 0;

        for (int j = 0; j < cool_dim; j++) {
            printf("pivot: %lf %p\n", pivot_addr->centroid[j], pivot_addr);
            printf("other: %lf %p\n", (*(triangulation + i*sizeof(simplex_t*)))->centroid[j], (*(triangulation + i*sizeof(simplex_t*))));
//            dist += pow(pivot_addr->centroid[j] - (*(triangulation + i*sizeof(simplex_t*)))->centroid[i], 2);
            dist += (pivot_addr->centroid[j] - (*(triangulation + i*sizeof(simplex_t*)))->centroid[j]) * (pivot_addr->centroid[j] - (*(triangulation + i*sizeof(simplex_t*)))->centroid[j]);
        }
        printf("%lf\n", dist);

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


