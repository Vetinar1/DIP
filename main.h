//
// Created by vetinari on 26.08.20.
//

#define cool_dim 2
#define qhull_dim 3

#define n_points 258
#define n_points 5

#ifndef MASTER_PROJECT_C_PART_MAIN_H
#define MASTER_PROJECT_C_PART_MAIN_H

#endif //MASTER_PROJECT_C_PART_MAIN_H

typedef struct tri tri_t;
struct tri {
    // Hierarchical triangle element representing a triangle from the Delaunay building process
    tri_t * children[3];        // Array of pointers to tris; Each triangle element has up to three children
    int points[cool_dim+1];     // Indices of the points of the tri
    int final;                  // Nonzero if this tri is part of the final triangulation

};

typedef struct qfacet qfacet_t;
struct qfacet {
    int points[qhull_dim];                  // points that make up the facet
    qfacet_t * neighbors[qhull_dim+1];      // Neighboring facets
    double normal[qhull_dim];               // Normal on the hyperplane
    int outside[n_points];                  // Need to keep this constant sized somehow so its a hash table
                                            // Contains all points "outside"/"above" this facet
                                            // TODO Should probably be a really long bit vector but i cant be bothered
    int outside_empty;                      // 1 if outside is filled with 0s, otherwise 0
    int in_visible, visited;                // Flags for algorithm
};

int contains(tri_t*, const double*, int);   // Checks if tri contains point
void invert_matrix(double*, int);
void calculate_normal(qfacet_t*, double*, int, int);
double signed_dist_to_facet(qfacet_t*, double*, int);