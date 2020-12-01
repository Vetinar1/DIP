//
// Created by vetinari on 26.11.20.
//


#ifndef MASTER_PROJECT_C_PART_COOLCOOL_H
#define MASTER_PROJECT_C_PART_COOLCOOL_H

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
#include "CoolPoint.h"
#include "CoolSimplex.h"
#include "CoolConst.h"


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
    Simplex<D> * btree;         // Points to the root of the simplex ball tree

    Simplex<D> * construct_simplex_btree_recursive(Simplex<D> **, int);
    Simplex<D> * find_nearest_neighbour_sbtree(Simplex<D> *, const double *, Simplex<D> *, double);

    int flips, interpolate_calls;
    int N_MAX, S_MAX;
public:
    Cool() {
        flips = 0;
        interpolate_calls = 0;
        avg_flips = 0;
        S_MAX = S;
        N_MAX = N;
        for (int i = 0; i < D; i++) {
            mins[i] = DBL_MAX;
            maxs[i] = -1 * DBL_MAX;
        }
    };
    double avg_flips;

    // Minimum and maximum values in each dimension
    double mins[D];
    double maxs[D];

    void reset();
    int read_files(std::string, std::string, std::string);
    void save_btree(std::string filename);
    int construct_btree();
    double interpolate(double * coords);
    void set_N_MAX(int n) { N_MAX = n; };
    void set_S_MAX(int s) { S_MAX = s; };
};


template<int N, int D, int S>
void Cool<N, D, S>::reset() {
    flips = 0;
    interpolate_calls = 0;
    avg_flips = 0;
    S_MAX = S;
    N_MAX = N;
    for (int i = 0; i < D; i++) {
        mins[i] = DBL_MAX;
        maxs[i] = -1 * DBL_MAX;
    }

    // Points[N] will be overwritten on next call to read_files()
}


template<int N, int D, int S>
int Cool<N, D, S>::construct_btree() {
    /**
     * This function serves as a public "adapter" to the actual ball tree construction function,
     * construct_point_btree_recursive().
     *
     * Returns 0 on success.
     */
    // Construct array of pointers
    Simplex<D> * simps[S];
    for (int i = 0; i < S_MAX; i++) {
        simps[i] = &(simplices[i]);
    }

    btree = construct_simplex_btree_recursive(simps, S_MAX);

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

        assert(dim_min <= dim_max);
        dim_spread = fabs(dim_max - dim_min);
        if (dim_spread > largest_spread) {
            largest_spread = dim_spread;
            lspread_dim    = i;
            avg            = (dim_max + dim_min) / 2;
        }
    }

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
void Cool<N, D, S>::save_btree(std::string filename) {
    /**
     * Saves the Ball tree. Format: [Index of simplex] [Index of left child] [Index of right child] [pbtree_radius_sq]
     *
     * Grows with O(N^2) ... TODO
     */
    std::ofstream file;
    file.open(filename);
    for (int i = 0; i < S_MAX; i++) {
        file << i << " ";
        int lchild = -1, rchild = -1;
        for (int j = 0; j < S_MAX; j++) {
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
double Cool<N, D, S>::interpolate(double * coords) {
    /**
     * Interpolate the given point.
     * First, find the closest simplex using the ball tree. Then, find the simplex containing the given point using
     * repeated "flips". Finally, interpolate using a weighted average (Delaunay).
     */

    Simplex<D> * best = NULL;
    Simplex<D> * nn = find_nearest_neighbour_sbtree(btree, coords, best, DBL_MAX);

    double * bary = nn->convert_to_bary(coords);
    int inside = nn->check_bary(bary);
    int flips = 0;

    while (!inside) {
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
        bary = nn->convert_to_bary(coords);
        inside = nn->check_bary(bary);

        flips++;
        if (flips > 100) {
            std::cerr << "Error: More than 100 flips." << std::endl;
            abort();
        }

    }

    // The actual interpolation step
    double val = 0;
    for (int i = 0; i < D+1; i++) {
        val += bary[i] * nn->points[i]->value;
    }

    delete[] bary;
    delete best;

    interpolate_calls++;
    flips += flips;
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


    /* Read points */
    file.open(cool_file);
    if (!file.is_open()) {
        std::cerr << "Error reading " << cool_file << std::endl;
        return 1;
    } else {
        std::cout << "Reading" << cool_file << std::endl;
    }

    int n = 0;
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
                mins[j] = points[i].coords[j];
            }
            if (points[i].coords[j] > maxs[j]) {
                maxs[j] = points[i].coords[j];
            }
        }

        n++;
        if (file.peek() == EOF) {
            break;
        }
    }
    N_MAX = n;


    file.close();

    /* Read Simplices */
    file.open(tri_file);
    if (!file.is_open()) {
        std::cerr << "Error reading " << tri_file << std::endl;
        return 2;
    } else {
        std::cout << "Reading" << tri_file << std::endl;
    }


    int s = 0;
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

        s++;
        if (file.peek() == EOF) {
            break;
        }
        if (skip == 1) {
            continue;
        }
        simplices[i].calculate_normals();
        simplices[i].validate_normals();
    }
    S_MAX = s;
    std::cout << "S: " << S << " S_MAX: " << S_MAX << std::endl;



    file.close();

    // Construct matrices for each simplex
    // https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Barycentric_coordinates_on_tetrahedra
    for (int s = 0; s < S_MAX; s++) {   // for each simplex

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
    } else {
        std::cout << "Reading" << neighbour_file << std::endl;
    }


    for (int i = 0; i < S_MAX; i++) {
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

    return 0;
}


#endif //MASTER_PROJECT_C_PART_COOLCOOL_H