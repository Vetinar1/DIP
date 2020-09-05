//
// Created by vetinari on 26.08.20.
//

#include "main.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

int main() {
    // T, nH, Values
    double points[n_points][cool_dim+1] = {
            {0, 1, 12},
            {1, 0, 34},
            {0, 0, 13},
            {0, 0.5, 34},
            {2, 3, 34},
    };
    tri_t tri1;
    tri1.points[0] = 0;
    tri1.points[1] = 1;
    tri1.points[2] = 2;

    int dummy = 0;
    dummy = contains(&tri1, &points[0][0], 3);
//    FILE * f = fopen("../data.csv", "r");
//    if (f == NULL) {
//        perror("failed: ");
//        return 1;
//    }
//    for (int i = 0; i < n_points; i++) {
//        fscanf(f, "%lf,%lf,%lf\n", &points[i][0], &points[i][1], &points[i][2]);
//    }
//
//    for (int i = 0; i < n_points; i++) {
//        printf("%f %f %f\n", points[i][0], points[i][1], points[i][2]);
//    }

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


int contains(tri_t * tri, const double * points, int pt_idx) {
    /*
     * Check if point in points at index pt_index is contained by tri_t tri.
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


void qhull(double * points, int N, int d) {
    /*
     * https://www.researchgate.net/publication/2641780_The_QuickHull_Algorithm_for_Convex_Hulls
     */
    int * unassigned_points = (int*)calloc(N, sizeof(int));
    for (int i = 0; i < N; i++) {
        unassigned_points[i] = 1;
    }

    int qft_s = sizeof(qfacet_t);   // qfacet_t size
    int facet_count = d + 1;        // start out with enough facets to build one simplex
    qfacet_t * hull = malloc(facet_count * qft_s);

    int * init_points = malloc((d+1)*sizeof(int));

    // 1. Create simplex of d+1 points
    // TODO select points with lowest/largest coordinates; for now just take the first few. DEGENERATE!
    // a) Select points
    for (int i = 0; i < d+1; i++) {
        init_points[i] = i;
    }

    // b) Build facets
    // TODO Build "constructor" function
    for (int i = 0; i < d+1; i++) {
        qfacet_t new_facet;
        // Each initial facet is made up of all initial points except one
        for (int j = 0; j < d; j++) {
            new_facet.points[j] = j < i ? init_points[j] : init_points[j+1];
        }
        // Initialize hash table of outside points to zeros
        for (int j = 0; j < N; j++) {
            new_facet.outside[j] = 0;
        }
        new_facet.outside_empty = 1;
        new_facet.visited = 0;
        new_facet.in_visible = 0;
        *(hull + i*qft_s) = new_facet;
    }

    // c) Link facets and calculate their normals
    // All initial facets are neighbors
    for (int i = 0; i < d+1; i++) {
        for (int j = 0; j < d+1; j++) {
            *(hull + i*qft_s)->neighbors[j] = j < i ? *(hull + j * qft_s) : *(hull + (j+1) * qft_s);
        }

        calculate_normal(hull + i*qft_s, points, N, d);
        // TODO: Verify normal points outwards
    }

    // 2. for each facet F
    for (int i = 0; i < d+1; i++) {
        // 2.1 for each unassigned point p
        for (int j = 0; j < N; j++) {
            // 2.2 if p is above F
            // "Above" if the signed distance is positive - dot product of plane normal with point;
            // equivalent to projection of point vector on normal vector
            double signed_dist = signed_dist_to_facet((hull + qft_s*i), (points + j*d), d);
            if (signed_dist > 0) {
                // 2.3 assign p to Fs outside set
                unassigned_points[j] = 0;
                (hull + qft_s*i)->outside[j] = 1;
                (hull + qft_s*i)->outside_empty = 0;
            }
        }
    }

    // 3. for each facet F with a non-empty outside set
    for (;;) {
        // Find a facet with a non empty outside set
        qfacet_t curr_facet;
        int curr_facet_exists = 0;
        for (int i = 0; i < facet_count; i++) {
            if ((hull + i*qft_s)->outside_empty == 0) {
                curr_facet = *(hull + i*qft_s);
                curr_facet_exists = 1;
                break;
            }
        }
        if (curr_facet_exists == 0) {
            break;
        }

        // 3.1 select the furthest point p of Fs outside set
        double largest_dist = 0;
        int largest_dist_index = -1;
        for (int i = 0; i < N; i++) {
            if (curr_facet.outside[i] == 1) {
                double signed_dist = signed_dist_to_facet(&curr_facet, (points + i*d), d);
                assert(signed_dist > 0);
                if (signed_dist > largest_dist) {
                    largest_dist = signed_dist;
                    largest_dist_index = i;
                }
            }
        }
        assert(largest_dist_index != -1);

        // 3.2 initialize the visible set V to F (array of pointers to structs)
        // The maximum number of visible simplices is all simplices minus 1
        qfacet_t ** visible = (qfacet_t **)malloc(facet_count * sizeof(qfacet_t *));
        int visible_count = 1;      // The actual number of visible simplices

        // The current facet is always visible
        visible[0] = &curr_facet;

        // 3.3 for all unvisited neighbors N of facets in V
        for (;;) {
            int new_visible_this_iteration = 0;
            // For all unvisited visible facets
            for (int i = 0; i < visible_count; i++) {
                if (visible[i]->visited == 1) {
                    continue;
                }
                visible[i]->visited = 1;

                // Determine for all neighbors not yet in "visible" if they are visible,
                // and if yes, add them to visible. They are unvisited
                for (int j = 0; j < qhull_dim+1; j++) {
                    if (visible[i]->neighbors[j]->in_visible == 1) {
                        continue;
                    }
                    // TODO: Can further improve this by keeping track of not-visible facets, but i *really* dont care rn
                    double signed_dist = signed_dist_to_facet(
                        visible[i]->neighbors[j],
                        (points + d*largest_dist_index),
                        d
                    );

                    // 3.3.1 if p is above N
                    if (signed_dist > 0) {
                        // 3.3.2 add N to V
                        visible[visible_count + new_visible_this_iteration + 1] = visible[i]->neighbors[j];
                        visible[i]->neighbors[j]->in_visible = 1;
                        new_visible_this_iteration++;
                    }
                }
            }

            if (new_visible_this_iteration == 0) {
                break;
            }
            visible_count += new_visible_this_iteration;
        }

        // 3.4 the boundary of V is the set of horizon ridges H
        // horizon ridge = ridge across which the neighbor is not visible
        qfacet_t * new_facets = (qfacet_t*)malloc(1000*qft_s);  // TODO how much memory do i need?
        for (int i = 0; i < visible_count; i++) {
            for (int j = 0; j < qhull_dim+1; j++) {
                if
            }
        }
    }

        // for each ridge R in H
            // create a new facet from R and p
            // link the new facet to its neighbors
        // for each new facet F'
            // for each unassigned point q in an outside set of a facet in V
                // if q is above F'
        // delete the facets in V
}


double signed_dist_to_facet(qfacet_t * qf, double * p, int d) {
    /*
     * calculate signed distance from facet struct to point found at pointer p. in space of dimension d
     */
    // Note that plane offset needs to be taken into account. (Draw a diagram)
    // Since its a plane, ANY point will work as offset, so i just use the first of the facet
    // TODO: Verify
    double * subtr;
    for (int k = 0; k < d; k++) {
        subtr[k] = p[k] - qf->points[0];
    }

    double signed_dist = 0;
    for (int k = 0; k < d; k++) {
        signed_dist += subtr[k] * qf->normal[k];
    }

    return signed_dist;
}


// TODO N and d mean different things here than in the other function, unify
void calculate_normal(qfacet_t * qf, double * points, int N, int d) {
    /*
     * Calculate normal of a qfacet
     * https://math.stackexchange.com/a/3398303
     *
     */
    // Build matrix
    double * matrix = (double*)calloc(N*N, sizeof(double));

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            *(matrix + N*i + j) = *(points + d*qf->points[i] + j);
        }

        qf->normal[i] = 1;
    }

    solve_linear(&matrix[0], &qf->normal[0], N);

    double length = 0;
    for (int i = 0; i < N; i++) {
        length += qf->normal[i] * qf->normal[i];
    }
    length = sqrt(length);
    for (int i = 0; i < N; i++) {
        qf->normal[i] /= length;
    }
}


