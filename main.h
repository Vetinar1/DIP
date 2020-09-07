//
// Created by vetinari on 26.08.20.
//

#define cool_dim 2

#define n_points 239
#define n_simplices 452

#ifndef MASTER_PROJECT_C_PART_MAIN_H
#define MASTER_PROJECT_C_PART_MAIN_H

#endif //MASTER_PROJECT_C_PART_MAIN_H

typedef struct simplex simplex_t;
struct simplex {
    // Hierarchical triangle element representing a triangle from the Delaunay building process
    int points[cool_dim+1];             // Indices of the points of the simplex
    simplex_t * neighbors[cool_dim+1];  // The ith neighbor is the neighbor opposite of the ith vertex

    // For use in ball tree algorithm
    double centroid[cool_dim];
    double btree_radius_sq;
    simplex_t * lchild;
    simplex_t * rchild;
    int part_of_tree;
};

/*typedef struct qfacet qfacet_t;
struct qfacet {
    int points[qhull_dim];                  // points that make up the facet
    qfacet_t * neighbors[qhull_dim+1];      // Neighboring facets
    double normal[qhull_dim];               // Normal on the hyperplane
    int outside[n_points];                  // Need to keep this constant sized somehow so its a hash table
                                            // Contains all points "outside"/"above" this facet
                                            // TODO Should probably be a really long bit vector but i cant be bothered
    int outside_empty;                      // 1 if outside is filled with 0s, otherwise 0
    int in_visible, visited;                // Flags for algorithm
};*/

int contains(simplex_t*, const double*, int);   // Checks if simplex contains point
void invert_matrix(double*, int);
//void calculate_normal(qfacet_t*, double*, int, int);
//double signed_dist_to_facet(qfacet_t*, double*, int);

simplex_t * build_ball_tree(int, simplex_t**);