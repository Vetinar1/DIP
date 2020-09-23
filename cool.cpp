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

template<int D> class Point;
template<int D> class Simplex;
template<int N, int D, int S> class Cool;

template<int D>
class Point {
    /**
     * Class representing a single data point. Mostly used to find simplices it is part of.
     */
private:
    double coords[D];
    double value;
    double btree_radius_sq;

public:
    std::vector<Simplex<D>> simplices;  // All simplices that contain this point
    Point<D> * lchild;
    Point<D> * rchild;

    Point() {
        btree_radius_sq = 0;
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
//    int neighbour_indices[D+1];         // One neighbour opposite every point
//    Simplex * neighbour_pointers[D+1];
    double centroid[D];
    double btree_radius_sq;

    void invert_T();
    double * convert_to_bary(const double *);
    int check_bary(const double *);
public:
    Simplex * lchild;
    Simplex * rchild;

    Simplex() {
        btree_radius_sq = 0;
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

    Point<D> * construct_point_btree_recursive(Point<D> **, int);
    Point<D> * find_nearest_neighbour_pbtree(Point<D> *, const double *, Point<D> *, double);
    Simplex<D> * construct_simplex_btree_recursive(Simplex<D> **, int);
    Simplex<D> * find_nearest_neighbour_sbtree(Simplex<D> *, const double *, Simplex<D> *, double);

    int flips, interpolate_calls;
public:
    Cool() {
//        points = new double[N][D+1];
//        simplices = new Simplex<D>[S];
        flips = 0;
        interpolate_calls = 0;
        avg_flips = 0;
    };
    double avg_flips;

    int read_files(std::string, std::string, std::string);
    void save_pbtree(std::string);
    void save_sbtree(std::string);
    int construct_point_btree();
    int construct_simplex_btree();
    double interpolate_sbtree(double *);
    double interpolate_pbtree(double *);
};


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
        pts[i] = &(simplices[i]);
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

    sbtree = construct_point_btree_recursive(simps, S);

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
        pts[0]->btree_radius_sq = 0;

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
    pivot_addr->btree_radius_sq = 0;
    for (int i = 0; i < n; i++) {
        double dist = 0;
        for (int j = 0; j < D; j++) {
            dist += pow(pts[i]->coords[j] - pivot_addr->coords[j], 2);
        }

        if (dist > pivot_addr->btree_radius_sq) {
            pivot_addr->btree_radius_sq = dist;
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
        simps[0]->btree_radius_sq = 0;

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
    pivot_addr->btree_radius_sq = 0;
    for (int i = 0; i < n; i++) {
        double dist = 0;
        for (int j = 0; j < D; j++) {
            dist += pow(simps[i]->centroid[j] - pivot_addr->centroid[j], 2);
        }

        if (dist > pivot_addr->btree_radius_sq) {
            pivot_addr->btree_radius_sq = dist;
        }
    }

    delete[] L;
    delete[] R;

    return pivot_addr;
}


template<int N, int D, int S>
void Cool<N, D, S>::save_pbtree(std::string filename) {
    /**
     * Saves the Ball tree. Format: [Index of simplex] [Index of left child] [Index of right child] [btree_radius_sq]
     *
     * Grows with O(N^2) ... TODO
     */
    std::ofstream file;
    file.open(filename);
    for (int i = 0; i < S; i++) {
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
        file << lchild << " " << rchild << " " << points[i].btree_radius_sq << std::endl;
    }
    file.close();
}


template<int N, int D, int S>
void Cool<N, D, S>::save_sbtree(std::string filename) {
    /**
     * Saves the Ball tree. Format: [Index of simplex] [Index of left child] [Index of right child] [btree_radius_sq]
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
        file << lchild << " " << rchild << " " << simplices[i].btree_radius_sq << std::endl;
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
    if (dist2 - root->btree_radius_sq >= min_dist2) {
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
    if (dist2 - root->btree_radius_sq >= min_dist2) {
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
    Simplex<D> * best = NULL;
    Simplex<D> * nn = find_nearest_neighbour_sbtree(sbtree, coords, best, DBL_MAX);

    double * bary = nn->convert_to_bary(coords);

    int inside = nn->check_bary(bary);

    double dist = 0;
    for (int i = 0; i < D; i++) {
        dist += pow(coords[i] - nn->centroid[i], 2);
    }
    dist = sqrt(dist);
    std::cout << "Dist: " << dist << std::endl;


    int dbg_count = 0;
    while (!inside) {
        // If the point is not contained in the simplex, the "most negative" barycentric coordinate denotes the one
        // "most opposite" of our coordinates. Take the simplex' neighbor on the opposite of that opposite,
        // and try again
        double min_bary = 1;
        int    min_bary_index = -1;
        for (int i = 0; i < D+1; i++) {
            if (bary[i] < min_bary) {
                min_bary = bary[i];
                min_bary_index = i;
            }
        }
        assert(min_bary < 1 && min_bary_index != -1);

        std::cout << "Flip " << dbg_count << std::endl;
        nn = nn->neighbour_pointers[min_bary_index];

        for (int i = 0; i < S; i++) {
            if (nn == &(simplices[i])) {
                std::cout << "Simplex: " << i << std::endl;
            }
        }

        double dist = 0;
        for (int i = 0; i < D; i++) {
            dist += pow(coords[i] - nn->centroid[i], 2);
        }
        dist = sqrt(dist);
        std::cout << "Dist: " << dist << std::endl;

        bary = nn->convert_to_bary(coords);

        inside = nn->check_bary(bary);

        dbg_count++;
        if (dbg_count > 100) {
            std::cerr << "Error: More than 100 flips." << std::endl;
            break;
        }
    }

    // The actual interpolation step
    double val = 0;
    for (int i = 0; i < D+1; i++) {
        val += bary[i] * nn->points[i][D];  // The Dth "coordinate" is the function value
    }

    delete[] bary;

    interpolate_calls++;
    flips += dbg_count;
    avg_flips = flips / (float) interpolate_calls;

    return val;
}


template<int N, int D, int S>
double Cool<N, D, S>::interpolate_pbtree(double * coords) {
    /**
     * Interpolate the given point.
     * Find closest point, then check all simplices touching that point
     */
    Point<D> * best = NULL;
    Point<D> * nn = find_nearest_neighbour(pbtree, coords, best, DBL_MAX);

    double * bary;
    int simp_index;

    for (int i = 0; i < nn->simplices.size()+1; i++) {
        if (i >= nn->simplices.size()) {
            std::cerr << "Error: Target not found in vicinity of nearest neighbor" << std::endl;
            abort();
        }

        int inside;
        bary = nn->simplices[i]->convert_to_bary(coords);
        inside = nn->simplices[i]->check_bary(bary);

        if (inside) {
            simp_index = i;
            break;
        }
    }

    // The actual interpolation step
    double val = 0;
    for (int i = 0; i < D+1; i++) {
        val += bary[i] * nn->simplices[simp_index]->points[i]->val;  // The Dth "coordinate" is the function value
    }

    delete[] bary;

    interpolate_calls++;
    flips += dbg_count;
    avg_flips = flips / (float) interpolate_calls;

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

    /* Read points.csv */
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
        points.simplices.shrink_to_fit();
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