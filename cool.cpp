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
    double coords[D];
    double value;
    double pbtree_radius_sq;

public:
    std::vector<Simplex<D>*> simplices;  // All simplices that contain this point
    Point<D> * lchild;
    Point<D> * rchild;

    Point() {
        pbtree_radius_sq = 0;
    }
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
    std::unordered_set<int> block_ids;         // Intersecting blocks

    void invert_T();
    double * convert_to_bary(const double *);
    int check_bary(const double *);
public:
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
//    double * points;
//    Simplex<D> * simplices;
    Simplex<D> * sbtree;         // Points to the root of the simplex ball tree
    Point<D> * pbtree;           // Points to the root of the point ball tree

    // D-dimensional array to pointers of vectors to pointers to vectors of pointers to int
    // Each dimension (D) has an arbitrary number (-> vector of pointers to vector) of pages (vector of int)
    // TODO: Somehow pass vector lengths as input
    std::vector<std::vector<int>*> * slices[D];
    int n_slices[D];
    int n_blocks[D];            // Number of block subdivisions in each direction. TODO: Make argument
    float block_size[D];        // Length of a block sized subdivision in dimension D. Used a lot -> Keep around

    double dims[D][2];          // Edges of dimensions; [min, max]

    Point<D> * construct_point_btree_recursive(Point<D> **, int);
    Point<D> * find_nearest_neighbour_pbtree(Point<D> *, const double *, Point<D> *, double);
    Simplex<D> * construct_simplex_btree_recursive(Simplex<D> **, int);
    Simplex<D> * find_nearest_neighbour_sbtree(Simplex<D> *, const double *, Simplex<D> *, double);
    int get_block_id(const double * );
    double * get_block_corners(int);
    void find_block_intersections(Simplex<D> *);

    int sbtree_flips, sbtree_interpolate_calls;
    int pbtree_flips, pbtree_interpolate_calls;
public:
    Cool() {
//        points = new double[N][D+1];
//        simplices = new Simplex<D>[S];
        sbtree_flips = 0;
        sbtree_interpolate_calls = 0;
        avg_sbtree_flips = 0;
        pbtree_flips = 0;
        pbtree_interpolate_calls = 0;
        for (int i = 0; i < D; i++) {
            n_blocks[i] = 100;
            n_slices[i] = 10;

            for (int j = 0; j < n_slices[i]; j++) {
                (*slices[i]).push_back(new std::vector<int>);
            }

            dims[i][0] = DBL_MAX;
            dims[i][1] = -1*DBL_MAX;
        }
    };
    double avg_sbtree_flips, avg_pbtree_flips;

    int read_files(std::string, std::string, std::string);
    void save_pbtree(std::string);
    void save_sbtree(std::string);
    int construct_point_btree();
    int construct_simplex_btree();
    double interpolate_sbtree(double *);
    double interpolate_pbtree(double *);
};


template<int N, int D, int S>
int Cool<N, D, S>::get_block_id(const double * coords) {
    // Assigns an id like
    // id = x + width * (y + depth * (z + height * (...)))
    int id = 0;
    for (int i = D-1; i >= 0; i--) {
        int ith_ind = (int) ((coords[i] - dims[i][0]) / block_size[i]);
        id = ith_ind + n_blocks[i] * id;
    }
}

template<int N, int D, int S>
double * Cool<N, D, S>::get_block_corners(int id) {
    // Retrieve corners of a block from id - inverse of get_block_id
    // Reconstruct from id = x + width * (y + depth * (z + height * (...)))
    // The corners are returned in no specific order
    double * corners = new double[D*pow(2, D)];     // 2^D points with D dimensions

    // Find first corner the one determined by the id
    for (int i = 0; i < D; i++) {
        corners[i] = id % n_blocks[i];
        id = (id - corners[i]) / n_blocks[i];
        corners[i] = dims[i][0] + corners[i] * block_size[i];
    }

    // Now find the other corners through repeated mirroring
    // TODO Verify this actually works
    int point_counter = 1;
    for (int i = 0; i < D; i++) {   // Mirror each dimension once
        for (int j = 0; i < point_counter; j++) {   // Mirror each point we have so far
            for (int k = 0; k < D; k++) {   // Each coordinate of those points
                if (k == i) {
                    corners[D*point_counter + D*j + k] = corners[D*j + k] + block_size[i];
                } else {
                    corners[D*point_counter + D*j + k] = corners[D*j + k];
                }
            }
        }
    }

    return corners;
}

template<int N, int D, int S>
void Cool<N, D, S>::find_block_intersections(Simplex<D> * s) {
    std::vector<int> discarded;
    std::vector<int> in_queue;

    s->block_ids.insert(get_block_id(s->points[0]->coords));
    in_queue.push_back(get_block_id(s->points[0]->coords));

    while (!in_queue.empty()) {
        int curr = in_queue.back();
        in_queue.pop_back();

        double * curr_corners = get_block_corners(curr);

        // For all faces of the simplex, check if any of the two corners are on opposite side of the plane
        // The check checks if two corners of the cube have different signs in the same barycentric coordinate
        // If yes, consider this simplex to be intersecting that hypercube
        // TODO CHECK IF THAT ACTUALLY WORKS
        // If that is the case, check the neighbors of that hypercube, else discard it
        double * corner_bary = new double[(D+1)];
        double * corner_bary_new = new double[(D+1)];
        for (int i = 0; i < pow(2, D); i++) {
            corner_bary_new = s->convert_to_bary(&(curr_corners[i*D]));

            if (i == 0) {
                for (int j = 0; j < D; j++) {
                    corner_bary[j] = corner_bary_new[j];
                }
                continue;
            }

            for (int j = 0; j < D; j++) {
                if (sgn(corner_bary[j]) != sgn(corner_bary_new[j])) {
                    s->block_ids.insert(curr);

                    // queue up more blocks
                    // TODO FLoating point inaccuracy might make this loop problematic
                    for (int k = 0; k < D; k++) {
                        int new1, new2;
                        double * new1_coords = new double[D];
                        double * new2_coords = new double[D];

                        for (int l = 0; l < D; l++) {
                            if (l != k) {
                                new1_coords[l] = curr_corners[l] + 0.5 * block_size[l];
                                new2_coords[l] = curr_corners[l] + 0.5 * block_size[l];
                            } else {
                                new1_coords[l] = curr_corners[l] + 1.5 * block_size[l];
                                new2_coords[l] = curr_corners[l] - 0.5 * block_size[l];
                            }
                        }

                        new1 = get_block_id(new1_coords);
                        new2 = get_block_id(new2_coords);

                        // I am so, so sorry
                        // https://stackoverflow.com/questions/3450860/check-if-a-stdvector-contains-a-certain-object
                        if (std::find(discarded.begin(), discarded.end(), new1) != discarded.end()) {
                        } else if (std::find(in_queue.begin(), in_queue.end(), new1) != in_queue.end()){
                        } else if (s->block_ids.count(new1) == 0) {
                            in_queue.push_back(new1);
                        }
                        if (std::find(discarded.begin(), discarded.end(), new2) != discarded.end()) {
                        } else if (std::find(in_queue.begin(), in_queue.end(), new2) != in_queue.end()){
                        } else if (s->block_ids.count(new2) == 0) {
                            in_queue.push_back(new2);
                        }


                        delete[] new1_coords;
                        delete[] new2_coords;
                    }
                    goto ENDOFWHILE;
                } else {
                    corner_bary[j] = corner_bary_new[j];
                }
            }
        }

        // This point is only reached if there is no intersection
        discarded.push_back(curr);

        ENDOFWHILE:
        delete[] corner_bary;
        delete[] curr_corners;
    }

    discarded.clear();
}


template<int N, int D, int S>
int Cool<N, D, S>::construct_point_btree() {
    /**
     * This function serves as a public "adapter" to the actual ball tree construction function,
     * construct_point_btree_recursive().
     *
     * Returns 0 on success.
     */
    // Construct array of pointers
    Point<D> * pts[N];
    for (int i = 0; i < N; i++) {
        pts[i] = &(points[i]);
    }

    pbtree = construct_point_btree_recursive(pts, N);

    return 0;
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
Point<D> * Cool<N, D, S>::construct_point_btree_recursive(Point<D> ** pts, int n) {
    /**
     * Simplex<D> pts       Pointer to array of pointers to Points to organize in tree
     * int n                Number of elements in array
     */

    // Input validation + what if theres only one element left?
    if (n < 0) {
        std::cerr << "Illegal argument in construct_point_btree: n = " << n << std::endl;
    } else if (n == 1) {
        pts[0]->lchild = NULL;
        pts[0]->rchild = NULL;
        pts[0]->pbtree_radius_sq = 0;

        return *pts;
    }

    // More than one element
    // 1. Find dimension of greatest spread
    double largest_spread = -1;
    int lspread_dim = -1;
    double avg = 0;             // Traditionally Ball trees work with the median, but that requires sorting TODO

    for (int i = 0; i < D; i++) {   // every dimension
        double dim_min =    DBL_MAX;
        double dim_max = -1*DBL_MAX;
        double dim_spread = 0;
        for (int j = 0; j < n; j++) {   // every point
            double val = pts[j]->coords[i];
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

    // 2. Find pivot - point closest to avg                          and
    // 3. Group points into sets to the left and right of average (L, R)
    double min_dist = DBL_MAX;
    Point<D> * pivot_addr = NULL;

    Point<D> ** L = new Point<D> * [n];
    Point<D> ** R = new Point<D> * [n];
    int lcount = 0;
    int rcount = 0;

    for (int i = 0; i < n; i++) {
        double val = pts[i]->coords[lspread_dim];
        Point<D> * addr = pts[i];
        double dist = fabs(avg - val);

        if (dist < min_dist) {
            // Add old pivot to L or R
            if (pivot_addr) {
                if (pivot_addr->coords[lspread_dim] <= avg) {
                    L[lcount] = pivot_addr;
                    lcount++;
                } else {
                    R[rcount] = pivot_addr;
                    rcount++;
                }
            }

            // Save new pivot
            min_dist = dist;
            pivot_addr = pts[i];

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
        pivot_addr->lchild = construct_point_btree_recursive(L, lcount);
        assert(pivot_addr->lchild != pivot_addr);
    } else {
        pivot_addr->lchild = NULL;
    }
    if (rcount > 0) {
        pivot_addr->rchild = construct_point_btree_recursive(R, rcount);
    } else {
        pivot_addr->rchild = NULL;
    }

    // This should only happen if there is only one element left, and in that case this line shouldnt be reached
    assert(pivot_addr->lchild != NULL || pivot_addr->rchild != NULL);

    // 5. Determine ball radius
    pivot_addr->pbtree_radius_sq = 0;
    for (int i = 0; i < n; i++) {
        double dist = 0;
        for (int j = 0; j < D; j++) {
            dist += pow(pts[i]->coords[j] - pivot_addr->coords[j], 2);
        }

        if (dist > pivot_addr->pbtree_radius_sq) {
            pivot_addr->pbtree_radius_sq = dist;
        }
    }

    delete[] L;
    delete[] R;

    return pivot_addr;
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

        assert(dim_min < dim_max);
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
void Cool<N, D, S>::save_pbtree(std::string filename) {
    /**
     * Saves the Ball tree. Format: [Index of simplex] [Index of left child] [Index of right child] [pbtree_radius_sq]
     *
     * Grows with O(N^2) ... TODO
     */
    std::ofstream file;
    file.open(filename);
    for (int i = 0; i < N; i++) {
        file << i << " ";
        int lchild = -1, rchild = -1;
        for (int j = 0; j < S; j++) {
            if (points[i].lchild == &points[j]) {
                lchild = j;
            }
            if (points[i].rchild == &points[j]) {
                rchild = j;
            }
        }
        file << lchild << " " << rchild << " " << points[i].pbtree_radius_sq << std::endl;
    }
    file.close();
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
Point<D> * Cool<N, D, S>::find_nearest_neighbour_pbtree(Point<D> * root, const double * target, Point<D> * best, double min_dist2) {
    /**
     * Recursive function to find the nearest neighbor of point target in tree root.
     * Adapted from https://en.wikipedia.org/wiki/Ball_tree#Pseudocode_2
     *
     * TODO I think this algorithm can be optimized to reduce the number of distance calculations
     *  -> talk to Klaus
     *
     * root             Root of subtree to search in
     * target           Coordinates of input point
     * best             Point with closest coords found so far
     * min_dist2        Distance to the coords of that simplex
     */
    double dist2 = 0;  // squared distance
    for (int i = 0; i < D; i++) {
        dist2 += pow(target[i] - root->coords[i], 2);
    }

    // Recursion exit condition - root is further than current closest neighbour
    if (dist2 - root->pbtree_radius_sq >= min_dist2) {
//        std::cout << "Worse." << std::endl;
        return best;
    }

    // assert(root->lchild != NULL || root->rchild != NULL);

    // Perhaps the current point is the new closest?
    if (dist2 < min_dist2) {
        min_dist2 = dist2;
        best = root;

//        int index = -1;
//        for (int i = 0; i < S; i++) {
//            if (best == &simplices[i]) {
//                index = i;
//                break;
//            }
//        }
//        std::cout << "New best: " << index << " dist: " << min_dist2 << std::endl;
    }

    // Find closest child
    double ldist2 = 0;
    double rdist2 = 0;

    if (root->lchild != NULL) {
        for (int i = 0; i < D; i++) {
            ldist2 += pow(target[i] - root->lchild->coords[i], 2);
        }
    }
    if (root->rchild != NULL) {
        for (int i = 0; i < D; i++) {
            rdist2 += pow(target[i] - root->rchild->coords[i], 2);
        }
    }

    // Recurse into closest child first
    // I wonder if this would be prettier with gotos...
    Point<D> * old_best = best;
    if (ldist2 <= rdist2) {
        if (root->lchild != NULL) {
//            std::cout << "Recursing left" << std::endl;
            best = find_nearest_neighbour_pbtree(root->lchild, target, best, min_dist2);
        }
        if (best != old_best) {
            // Recalculate min_dist2
            // TODO: Inefficient.
            min_dist2 = 0;
            for (int i = 0; i < D; i++) {
                min_dist2 += pow(target[i] - best->coords[i], 2);
            }
//            std::cout << "New min_dist2: " << min_dist2 << std::endl;
        }
        if (root->rchild != NULL) {
//            std::cout << "Recursing right" << std::endl;
            best = find_nearest_neighbour_pbtree(root->rchild, target, best, min_dist2);
        }
    } else {
        if (root->rchild != NULL) {
//            std::cout << "Recursing left" << std::endl;
            best = find_nearest_neighbour_pbtree(root->rchild, target, best, min_dist2);
        }
        if (best != old_best) {
            // Recalculate min_dist2
            // TODO: Inefficient.
            min_dist2 = 0;
            for (int i = 0; i < D; i++) {
                min_dist2 += pow(target[i] - best->coords[i], 2);
            }
//            std::cout << "New min_dist2: " << min_dist2 << std::endl;
        }
        if (root->lchild != NULL) {
//            std::cout << "Recursing right" << std::endl;
            best = find_nearest_neighbour_pbtree(root->lchild, target, best, min_dist2);
        }
    }
//    std::cout << "Going up " << std::endl;
//    std::cout << std::endl;

    return best;
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
//        std::cout << "Worse." << std::endl;
        return best;
    }

    // assert(root->lchild != NULL || root->rchild != NULL);

    // Perhaps the current point is the new closest?
    if (dist2 < min_dist2) {
        min_dist2 = dist2;
        best = root;

//        int index = -1;
//        for (int i = 0; i < S; i++) {
//            if (best == &simplices[i]) {
//                index = i;
//                break;
//            }
//        }
//        std::cout << "New best: " << index << " dist: " << min_dist2 << std::endl;
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
//            std::cout << "Recursing left" << std::endl;
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
//            std::cout << "Recursing right" << std::endl;
            best = find_nearest_neighbour_sbtree(root->rchild, target, best, min_dist2);
        }
    } else {
        if (root->rchild != NULL) {
//            std::cout << "Recursing left" << std::endl;
            best = find_nearest_neighbour_sbtree(root->rchild, target, best, min_dist2);
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
        if (root->lchild != NULL) {
//            std::cout << "Recursing right" << std::endl;
            best = find_nearest_neighbour_sbtree(root->lchild, target, best, min_dist2);
        }
    }
//    std::cout << "Going up " << std::endl;
//    std::cout << std::endl;

    return best;
}


template<int N, int D, int S>
double Cool<N, D, S>::interpolate_sbtree(double * coords) {
    /**
     * Interpolate the given point.
     * First, find the closest simplex using the ball tree. Then, find the simplex containing the given point using
     * repeated "flips". Finally, interpolate using a weighted average (Delaunay).
     */
    std::map<Simplex<D> *, int> visited;
    Simplex<D> * best = NULL;
    Simplex<D> * nn = find_nearest_neighbour_sbtree(sbtree, coords, best, DBL_MAX);

    double * bary = nn->convert_to_bary(coords);
    int inside = nn->check_bary(bary);

    double dist = 0;
    for (int i = 0; i < D; i++) {
        dist += pow(coords[i] - nn->centroid[i], 2);
    }
//    dist = sqrt(dist);
//    std::cout << "Coordinates: " << coords[0] << " " << coords[1] << std::endl;
//    for (int i = 0; i < S; i++) {
//        if (nn == &(simplices[i])) {
//            std::cout << "Simplex: " << i << " " << simplices[i].centroid[0] << " " << simplices[i].centroid[1] << std::endl;
//        }
//    }
//    std::cout << "Dist: " << dist << std::endl;


    int dbg_count = 0;
    while (!inside) {
        std::cout << dbg_count << std::endl;
        std::cout << "Coords: " << coords[0] << " " << coords[1] << " " << coords[2] << " " << coords[3] << std::endl;
        std::cout << "Centroid: " << nn->centroid[0] << " " << nn->centroid[1] << " " << nn->centroid[2] << " " << nn->centroid[3] << std::endl;
        double dist = 0;
        for (int i = 0; i < D; i++) {
            dist += pow(nn->centroid[i] - coords[i], 2);
        }
        dist = sqrt(dist);
        std::cout << "Dist: " << dist << std::endl;
        // If the point is not contained in the simplex, the "most negative" barycentric coordinate denotes the one
        // "most opposite" of our coordinates. Take the simplex' neighbor on the opposite of that opposite,
        // and try again
        // If we have visited once before, check the second furthest coordinate instead etc.
        // Note: This could technically lead to a case where we have walk from an N-1 times visited simplex to
        // an N times visited simplex and would have to instead "go back". I think this case is negligible
//        int n_visits = 0;
//        if (visited.find(nn) != visited.end()) {
//            n_visits = visited[nn];
//        } else {
//            visited[nn] = 0;
//        }
//        if (n_visits >= D+1) {
//            std::cerr << "Error: Simplex has been visited maximum number of times before" << std::endl;
//        }

//        std::cout << "n_visits: " << n_visits << std::endl;
        // TODO: Better variable names
        double min_bary = -1 * DBL_MAX;
        int    min_bary_index = -1;
//        for (int i = 0; i < n_visits+1; i++) {
        for (int i = 0; i < 1; i++) {
//            std::cout << i << std::endl;
            double nth_min_bary = 1;
            double nth_min_bary_index = -1;
            for (int j = 0; j < D+1; j++) {
                if (bary[j] < nth_min_bary && bary[j] > min_bary) {
                    nth_min_bary = bary[j];
                    nth_min_bary_index = j;
                }
            }
            min_bary = nth_min_bary;
            min_bary_index = nth_min_bary_index;
        }
        assert(min_bary > (-1*DBL_MAX) && min_bary_index != -1);

//        std::cout << "Flip " << dbg_count << std::endl;

        visited[nn] += 1;
        nn = nn->neighbour_pointers[min_bary_index];

//        for (int i = 0; i < S; i++) {
//            if (nn == &(simplices[i])) {
//                std::cout << "Simplex: " << i << " " << simplices[i].centroid[0] << " " << simplices[i].centroid[1] << std::endl;
//            }
//        }

//        double dist = 0;
//        for (int i = 0; i < D; i++) {
//            dist += pow(coords[i] - nn->centroid[i], 2);
//        }
//        dist = sqrt(dist);
//        std::cout << "Dist: " << dist << std::endl;

        bary = nn->convert_to_bary(coords);

        inside = nn->check_bary(bary);

        dbg_count++;
        if (dbg_count > 100) {
            std::cerr << "Error: More than 100 flips." << std::endl;
            exit(1);
            break;
        }
    }

    // The actual interpolation step
    double val = 0;
    for (int i = 0; i < D+1; i++) {
        val += bary[i] * nn->points[i]->coords[D];  // The Dth "coordinate" is the function value
    }

    delete[] bary;

    sbtree_interpolate_calls++;
    sbtree_flips += dbg_count;
    avg_sbtree_flips = sbtree_flips / (float) sbtree_interpolate_calls;

    return val;
}


template<int N, int D, int S>
double Cool<N, D, S>::interpolate_pbtree(double * coords) {
    /**
     * Interpolate the given point.
     * Find closest point, then check all simplices touching that point
     */
    Point<D> * best = NULL;
    Point<D> * nnp = find_nearest_neighbour_pbtree(pbtree, coords, best, DBL_MAX);
    std::cout << "Simplices: " << nnp->simplices[0] << std::endl;

    Simplex<D> * nns = NULL;
    double min_dist2 = DBL_MAX;
    for (int i = 0; i < nnp->simplices.size(); i++) {
        std::cout << i << std::endl;
        double dist2 = 0;
        for (int j = 0; j < D; j++) {
            dist2 += pow(nnp->simplices[i]->centroid[j] - coords[j], 2);
        }
        if (dist2 < min_dist2) {
            min_dist2 = dist2;
            nns = nnp->simplices[i];
        }
    }
    assert(min_dist2 != DBL_MAX);

    std::map<Simplex<D> *, int> visited;


    double * bary = nns->convert_to_bary(coords);
    int inside = nns->check_bary(bary);

    double dist = 0;
    for (int i = 0; i < D; i++) {
        dist += pow(coords[i] - nns->centroid[i], 2);
    }
//    dist = sqrt(dist);
//    std::cout << "Coordinates: " << coords[0] << " " << coords[1] << std::endl;
//    for (int i = 0; i < S; i++) {
//        if (nn == &(simplices[i])) {
//            std::cout << "Simplex: " << i << " " << simplices[i].centroid[0] << " " << simplices[i].centroid[1] << std::endl;
//        }
//    }
//    std::cout << "Dist: " << dist << std::endl;


    int dbg_count = 0;
    while (!inside) {
//        std::cout << "Coords: " << coords[0] << " " << coords[1] << " " << coords[2] << " " << coords[3] << std::endl;
//        std::cout << "Centroid: " << nns->centroid[0] << " " << nns->centroid[1] << " " << nns->centroid[2] << " " << nns->centroid[3] << std::endl;
//        double dist = 0;
//        for (int i = 0; i < D; i++) {
//            dist += pow(nns->centroid[i] - coords[i], 2);
//        }
//        dist = sqrt(dist);
//        std::cout << "Dist: " << dist << std::endl;
        // If the point is not contained in the simplex, the "most negative" barycentric coordinate denotes the one
        // "most opposite" of our coordinates. Take the simplex' neighbor on the opposite of that opposite,
        // and try again
        // If we have visited once before, check the second furthest coordinate instead etc.
        // Note: This could technically lead to a case where we have walk from an N-1 times visited simplex to
        // an N times visited simplex and would have to instead "go back". I think this case is negligible
        int n_visits = 0;
        if (visited.find(nns) != visited.end()) {
            n_visits = visited[nns];
        } else {
            visited[nns] = 0;
        }
        if (n_visits >= D+1) {
            std::cerr << "Error: Simplex has been visited maximum number of times before" << std::endl;
            break;
        }

//        std::cout << "n_visits: " << n_visits << std::endl;
        // TODO: Better variable names
        double min_bary = -1 * DBL_MAX;
        int    min_bary_index = -1;
        std::cout << "dbg_count: " << dbg_count << std::endl;
        std::cout << "n_visits: " << n_visits << std::endl;
        for (int i = 0; i < n_visits+1; i++) {
//            std::cout << i << std::endl;
            double nth_min_bary = DBL_MAX;
            int nth_min_bary_index = -1;
            for (int j = 0; j < D+1; j++) {
                std::cout << bary[j] << std::endl;
                if (bary[j] < nth_min_bary && bary[j] > min_bary) {
                    nth_min_bary = bary[j];
                    nth_min_bary_index = j;
                }
            }
            min_bary = nth_min_bary;
            min_bary_index = nth_min_bary_index;
            std::cout << nth_min_bary << " " << min_bary << std::endl;
        }
        std::cout << std::endl;
        assert(min_bary > (-1*DBL_MAX) && min_bary_index != -1);
//        std::cout << std::endl;

//        std::cout << "Flip " << dbg_count << std::endl;

        visited[nns] += 1;
        nns = nns->neighbour_pointers[min_bary_index];

//        for (int i = 0; i < S; i++) {
//            if (nn == &(simplices[i])) {
//                std::cout << "Simplex: " << i << " " << simplices[i].centroid[0] << " " << simplices[i].centroid[1] << std::endl;
//            }
//        }

//        double dist = 0;
//        for (int i = 0; i < D; i++) {
//            dist += pow(coords[i] - nn->centroid[i], 2);
//        }
//        dist = sqrt(dist);
//        std::cout << "Dist: " << dist << std::endl;

        bary = nns->convert_to_bary(coords);

        inside = nns->check_bary(bary);

        dbg_count++;
        if (dbg_count > 100) {
            std::cerr << "Error: More than 100 flips." << std::endl;
            break;
        }
    }

    // The actual interpolation step
    double val = 0;
    for (int i = 0; i < D+1; i++) {
        val += bary[i] * nns->points[i]->coords[D];  // The Dth "coordinate" is the function value
    }

    delete[] bary;

    pbtree_interpolate_calls++;
    pbtree_flips += dbg_count;
    avg_pbtree_flips = pbtree_flips / (float) pbtree_interpolate_calls;

    return val;

/*    double * bary;
    int simp_index;

//    std::cout << "Interpolation coords: " << coords[0] << " " << coords[1] << std::endl;
//    int nni;
//    for (int i = 0; i < N; i++) {
//        if (&points[i] == nn) {
//            nni = i;
//            break;
//        }
//    }
//    std::cout << "NN: " << nni << " " << nn->coords[0] << " " << nn->coords[1] << std::endl;

    int dbg_count = 0;
    for (int i = 0; i < nn->simplices.size()+1; i++) {
        for (int j = 0; j < S; j++) {
            if (nn->simplices[i] == &simplices[j]) {
                std::cout << "Simplex " << i << ": " << j << " " << simplices[j].centroid[0] << " " << simplices[j].centroid[1] << std::endl;
                break;
            }
        }

        if (i >= nn->simplices.size()) {
            std::cerr << "Error: Target not found in vicinity of nearest neighbor" << std::endl;
            abort();
        }

        int inside;
        bary = nn->simplices[i]->convert_to_bary(coords);
        inside = nn->simplices[i]->check_bary(bary);

        dbg_count++;
        if (inside) {
            simp_index = i;
            std::cout << "Interpolation success" << std::endl << std::endl;
            break;
        }
    }

    // The actual interpolation step
    double val = 0;
    for (int i = 0; i < D+1; i++) {
        val += bary[i] * nn->simplices[simp_index]->points[i]->value;  // The Dth "coordinate" is the function value
    }

    delete[] bary;

    pbtree_interpolate_calls++;
    pbtree_tries += dbg_count;
    avg_pbtree_tries = pbtree_tries / (float) pbtree_interpolate_calls;

    return val;*/
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
        for (int j = 0; j < D+1; j++) {     // D coordinates, 1 value
            std::getline(linestream, value, ',');
            points[i].coords[j] = std::stod(value);

            if (points[i].coords[j] < dims[j][0]) {
                dims[j][0] = points[i].coords[j];
            }
            if (points[i].coords[j] > dims[j][1]) {
                dims[j][1] = points[i].coords[j];
            }
        }
    }

    for (int i = 0; i < D; i++) {
        assert(dims[i][1] > dims[i][0]);
        block_size[i] = (dims[i][1] - dims[i][0]) / n_blocks[i];
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
            points[std::stoi(value)].simplices.push_back(&simplices[i]);
        }

        for (int j = 0; j < D+1; j++) {
            for (int k = 0; k < D+1; k++) {
                if (j != k && buffer[j] == buffer[k]) {
                    std::cerr << "Warning: Degenerate triangle (index " << i << ", points " << buffer[j] << " and " << buffer[k] << ")" << std::endl;
                }
            }
        }
        simplices[i].calculate_centroid();

        // Assign simplex to slices
        for (int j = 0; j < D; j++) {           // for each dimension
            // find min and max of vertex coordinates in that dimension
            double d_min = DBL_MAX;
            double d_max = 0;
            for (int k = 0; k < D+1; k++) {     // for each vertex
                if (simplices[i].points[k].coords[j] < d_min) {
                    d_min = simplices[i].points[k].coords[j];
                }
                if (simplices[i].points[k].coords[j] > d_max) {
                    d_max = simplices[i].points[k].coords[j];
                }
            }

            assert(d_max > d_min);

            // Assign to respective slices
            for (int k = 0; k < n_slices[j]; k++) {
                // slice_val = dimension min + k * slice_width
                // slice_width = (dimension max - dimension min) / n_slices[j]
                double slice_val = dims[j][0] + k * (d_max - d_min) / n_slices[j];

                if (d_min < slice_val) {
                    continue;
                } else if (d_max > slice_val) {
                    continue;
                } else {
                    // Insert index of simplex - slices are automatically sorted!
                    (*slices[i])[k]->push_back(i);
                }
            }
        }

        // Assign blocks to simplex
        // Strategy: Take a block that is known to intersect with the simplex. Check all its neighboring blocks.
        // For the ones that also intersect, recurse until no more intersecting blocks are left.
        find_block_intersections(simplices[i]);
    }

    // Slices will not be modified from here on
    for (int i = 0; i < D; i++) {
        slices[i]->shrink_to_fit();
        for (int j = 0; j < slices[i]->size(); j++) {
            slices[i]->at(j)->shrink_to_fit();
        }
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

    // Shrink point vectors
    for (int i = 0; i < N; i++) {
        points[i].simplices.shrink_to_fit();
    }
    return 0;
}


template<int D>
void Simplex<D>::invert_T() {
    /**
     * Invert matrix in place using Gauss Jordan Elimination with partial pivoting
     * Adapated from Numerical Recipes, 2.1
     *
     * TODO Efficiency? This should only be called once per simplex, probably not an issue
     * TODO Adopt implicit pivoting since we are not in a normalized space
     *
     * TODO This was mostly copied from my older C version of this code, and has not been properly tested for C++
     */

    int row_idx [D];
    for (int i = 0; i < D; i++) {
        row_idx[i] = i;
    }

    // yes, you can work in place instead, but i cant be bothered right now
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


    for (int i = 0; i < D; i++) {
        // find pivot
        int pivot_row_idx = i;
        double pivot = T_inv[row_idx[i]][i];
        for (int j = i; j < D; j++) {
            if (fabs(T_inv[row_idx[j]][i]) > fabs(pivot)) {
                pivot = T_inv[row_idx[j]][i];
                pivot_row_idx = j;
            }
        }

        // pivoting
        int temp = row_idx[i];
        row_idx[i] = row_idx[pivot_row_idx];
        row_idx[pivot_row_idx] = temp;

        // Normalization of current row
        double norm_factor = T_inv[row_idx[i]][i];
        for (int j = 0; j < D; j++) {
            T_inv[row_idx[i]][j] /= norm_factor;
            unit[row_idx[i]][j] /= norm_factor;
        }

        // Subtracting the current row from all other rows
        for (int j = 0; j < D; j++) {
            if (j != i) {
                double factor = T_inv[row_idx[j]][i];
                for (int k = 0; k < D; k++) {
                    T_inv[row_idx[j]][k] -= factor * T_inv[row_idx[i]][k];
                    unit[row_idx[j]][k] -= factor * unit[row_idx[i]][k];
                }
            }
        }
    }

    for (int i = 0; i < D; i++) {
        for (int j = 0; j < D; j++) {
            T_inv[row_idx[i]][j] = unit[i][j]; // TODO verify
        }
    }
}


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
        if (bary[i] <= 0) {
            // Consider points.csv on edges to not be contained, in order to avoid degeneracies!
            return 0;
        }
        if (isnanf(bary[i])) {
            std::cerr << "Error: Barycentric coordinate " << i << " is nan" << std::endl;
            return 0;
        }
        bsum += bary[i];
    }

    if (bsum > 1) {
        return 0;
    }

    return 1;
}