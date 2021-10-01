//
// Created by vetinari on 14.12.20.
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
#include <unordered_map>
#include <algorithm>
#include <bitset>
#include "CoolPoint.h"
#include "CoolSimplex.h"
#include "CoolConst.h"
#include "CoolCool.h"

void Cool::reset() {
    /**
     * Resets all properties of the cool object (and associated simplex objects) so it can be reused. Intended
     * to be a preparation for read_files()
     */
    flips = 0;
    interpolate_calls = 0;
    avg_flips = 0;
    nullpointers_encountered = 0;
    for (int i = 0; i < D; i++) {
        mins[i] = DBL_MAX;
        maxs[i] = -1 * DBL_MAX;
    }

    // Delete ball tree
    // Note: All other simplex properties will be overwritten on next call to read_files()
    for (int s = 0; s < S_MAX; s++) {
        simplices[s].sbtree_radius_sq = 0;
        simplices[s].lchild = nullptr;
        simplices[s].rchild = nullptr;
    }
    btree = nullptr;

    S_LIM = S_MAX;
    N_LIM = DIP_NMAX;

    // Points[N] will be overwritten on next call to read_files()
}


int Cool::construct_btree() {
    /**
     * This function serves as a public "adapter" to the actual ball tree construction function,
     * construct_simplex_btree_recursive().
     *
     * Returns 0 on success.
     */
    // Construct array of pointers; allocate on heap to avoid stack overflows
    Simplex ** simps = new Simplex *[S_LIM];
    for (int i = 0; i < S_LIM; i++) {
        simps[i] = &(simplices[i]);
    }

    btree = construct_simplex_btree_recursive(simps, S_LIM);

    delete[] simps;

    return 0;
}


Simplex * Cool::construct_simplex_btree_recursive(Simplex ** simps, int n) {
    /**
     * Implements ball tree construction algorithm. For a given array of simplices, determine the "middle" one (pivot).
     * Split the array in half at this simplex, and call function recursively on each half. Determine ball tree radius
     * of current pivot, and hook it up to the returns of the left/right half function calls. Return pointer to current
     * pivot.
     *
     * Simplex simps     Pointer to array of pointers to simplices to organize in tree
     * int n             Number of elements in array
     */

    // Input validation + what if theres only one element left?
    if (n < 0) {
        std::cerr << "Illegal argument in construct_btree: n = " << n << std::endl;
    } else if (n == 1) {
        simps[0]->lchild = nullptr;
        simps[0]->rchild = nullptr;
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

#ifdef DIP_BALLTREE_CHECKS
        assert(dim_min <= dim_max);
#endif
        dim_spread = fabs(dim_max - dim_min);
        if (dim_spread > largest_spread) {
            largest_spread = dim_spread;
            lspread_dim    = i;
            avg            = (dim_max + dim_min) / 2;
        }
    }

#ifdef DIP_BALLTREE_CHECKS
    assert(largest_spread != -1 && lspread_dim != -1);
#endif

    // 2. Find pivot - simplex with centroid closest to avg          and
    // 3. Group points into sets to the left and right of average (L, R)
    double min_dist = DBL_MAX;
    Simplex * pivot_addr = nullptr;

    Simplex ** L = new Simplex * [n];
    Simplex ** R = new Simplex * [n];
    int lcount = 0;
    int rcount = 0;

    for (int i = 0; i < n; i++) {
        double val = simps[i]->centroid[lspread_dim];
        Simplex * addr = simps[i];
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

#ifdef DIP_BALLTREE_CHECKS
    assert(pivot_addr != nullptr);
    assert(rcount != 0 || lcount != 0);
    assert(min_dist < DBL_MAX);
    for (int i = 0; i < lcount; i++) {
        assert(L[i] != pivot_addr);
    }
    for (int i = 0; i < rcount; i++) {
        assert(R[i] != pivot_addr);
    }
#endif

    // 4. Recurse on L and R
    if (lcount > 0) {
        pivot_addr->lchild = construct_simplex_btree_recursive(L, lcount);
    } else {
        pivot_addr->lchild = NULL;
    }
    if (rcount > 0) {
        pivot_addr->rchild = construct_simplex_btree_recursive(R, rcount);
    } else {
        pivot_addr->rchild = NULL;
    }

#ifdef DIP_BALLTREE_CHECKS
    // This should only happen if there is only one element left, and in that case this line shouldnt be reached
    assert(pivot_addr->lchild != NULL || pivot_addr->rchild != NULL);
#endif

    // 5. Determine ball radius
    // squared distance to farthest element in subtree
    pivot_addr->sbtree_radius_sq = 0;
    for (int i = 0; i < n; i++) {   // n elements of input array
        double dist = 0;
        // distance calculation
        for (int j = 0; j < D; j++) {   // D coordinates
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


void Cool::save_btree(std::string filename) {
    /**
     * Saves the Ball tree. Format: [Index of simplex] [Index of left child] [Index of right child] [pbtree_radius_sq]
     *
     * Warning, this function was only used for debug purposes and may not work fast or correctly.
     *
     * std::string filename     Filename.
     */
    // Map simplex pointer -> index
    std::unordered_map<Simplex *, int> s_indices;
    for (int i = 0; i < S_LIM; i++) {
        s_indices[&(simplices[i])] = i;
    }
    std::ofstream file;
    file.open(filename);
    for (int i = 0; i < S_LIM; i++) {
        file << i << s_indices[simplices[i].lchild] << " "
             << s_indices[simplices[i].rchild] << " " << simplices[i].sbtree_radius_sq << std::endl;
    }
    file.close();
}


Simplex * Cool::find_nearest_neighbor_sbtree(Simplex * root, const double * target, Simplex * best, double min_dist2) {
    /**
     * Recursive function to find the nearest neighbor of point target in tree root.
     * Adapted from https://en.wikipedia.org/wiki/Ball_tree#Pseudocode_2
     *
     * root             Root of subtree to search in
     * target           Coordinates of input point
     * best             Simplex with closest centroid found so far
     * min_dist2        Distance to the centroid of that simplex
     */
    // Distance to current node
    double dist2 = 0;
    for (int i = 0; i < D; i++) {
        dist2 += pow(target[i] - root->centroid[i], 2);
    }

    // Recursion exit condition - target point is further outside the current node's ball than the distance
    // to the current closest neighbor -> there can't be a closer neighbor in this ball/subtree
    if (sqrt(dist2) - sqrt(root->sbtree_radius_sq) >= sqrt(min_dist2)) {
        return best;
    }

    // Perhaps the current point is the new closest?
    if (dist2 < min_dist2) {
        min_dist2 = dist2;
        best = root;
    }

    // Find closest child
    double ldist2 = 0;
    double rdist2 = 0;

    if (root->lchild != nullptr) {
        for (int i = 0; i < D; i++) {
            ldist2 += pow(target[i] - root->lchild->centroid[i], 2);
        }
    }
    if (root->rchild != nullptr) {
        for (int i = 0; i < D; i++) {
            rdist2 += pow(target[i] - root->rchild->centroid[i], 2);
        }
    }

    // Recurse into closest child first
    // I wonder if this would be prettier with gotos...
    Simplex * old_best = best;
    if (ldist2 <= rdist2) {
        if (root->lchild != nullptr) {
            best = find_nearest_neighbor_sbtree(root->lchild, target, best, min_dist2);
        }
        if (best != old_best) {
            // Recalculate min_dist2
            // TODO: Inefficient.
            min_dist2 = 0;
            for (int i = 0; i < D; i++) {
                min_dist2 += pow(target[i] - best->centroid[i], 2);
            }
        }
        if (root->rchild != nullptr) {
            best = find_nearest_neighbor_sbtree(root->rchild, target, best, min_dist2);
        }
    } else {
        if (root->rchild != nullptr) {
            best = find_nearest_neighbor_sbtree(root->rchild, target, best, min_dist2);
        }
        if (best != old_best) {
            // Recalculate min_dist2
            // TODO: Inefficient.
            min_dist2 = 0;
            for (int i = 0; i < D; i++) {
                min_dist2 += pow(target[i] - best->centroid[i], 2);
            }
        }
        if (root->lchild != nullptr) {
            best = find_nearest_neighbor_sbtree(root->lchild, target, best, min_dist2);
        }
    }

    return best;
}


double * Cool::interpolate(double * coords) {
    /**
     * Interpolate the given point.
     * First, find the closest simplex using the ball tree. Then, find the simplex containing the given point using
     * repeated "flips". Finally, interpolate using a weighted average (Delaunay).
     *
     * double * coords      Pointer to array containing the target coordinates. Length needs to match dimensionality D.
     *                      Contents need to match the files that have been read. Do not attempt to interpolate outside
     *                      of sampled area of parameter space.
     * returns              Pointer to array of interpolated values, length DIP_VARNR
     */

    // Warning! Input coordinates are modified
    for (int i = 0; i < D; i++) {
        if (coords[i] > CLAMP_MAX[i]) {
            coords[i] = CLAMP_MAX[i];
        } else if (coords[i] < CLAMP_MIN[i]) {
            coords[i] = CLAMP_MIN[i];
        }
    }

    // 1. Find closest simplex (by centroid) using the ball tree
    Simplex * best = nullptr;
    Simplex * nn = find_nearest_neighbor_sbtree(btree, coords, best, DBL_MAX);
    if (nn == nullptr) {
        std::cerr << "DIP Important WARNING: Nullpointer returned from ball tree!" << std::endl;
//        std::cerr << "This should NEVER happen, are you SURE everything is working fine?" << std::endl;
        return nullptr;
    }
    assert(nn != nullptr);

    // 2. Find simplex actually containing the target point using simplex flipping algorithm
    // Check if current nearest neighbor contains point
    double * bary = nn->convert_to_bary(coords);
    int inside = nn->check_bary(bary);
    int n_flips = 0;

    while (!inside and n_flips < MAX_FLIPS) {
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

            double dot = 0;
            for (int j = 0; j < D; j++) {
                dot += diff_vec[j] * nn->normals[i][j] / diff_len;
            }

            if (dot > best_dir_dot) {
                best_dir_dot = dot;
                best_dir = i;
            }
        }

        if (nn->neighbor_pointers[best_dir] == nullptr) {
            nullpointers_encountered++;
            std::cerr << "Warning: Encountered nullpointer in simplex traversal (flip " << n_flips << ")" << std::endl;
            std::cerr << "So far this happened " << nullpointers_encountered << " times" << std::endl;
            std::cerr << "Coordinates: ";
            for (int i = 0; i < D; i++) {
                std::cerr << coords[i] << " ";
            }
            std::cerr << std::endl;
            std::cerr << "Bary Coordinates: ";
            for (int i = 0; i < D; i++) {
                std::cerr << bary[i] << " ";
            }
            std::cerr << std::endl;
            std::cerr << "Simplex at address " << nn << ", neighbor " << best_dir << std::endl;
            nn->print_error_info();
            break;
        }
        // Check if current nearest neighbor contains target point
        nn = nn->neighbor_pointers[best_dir];
        delete[] bary;
        bary = nn->convert_to_bary(coords);
        inside = nn->check_bary(bary);

        n_flips++;
    }

    // The actual interpolation step
    double * interp_vals = new double[DIP_VARNR];
    for (int i = 0; i < DIP_VARNR; i++) {
        interp_vals[i] = 0;
        for (int j = 0; j < DIP_DIMS+1; j++) {
            interp_vals[i] += bary[i] * nn->points[j]->value[i];
        }
    }
    
    delete[] bary;

#ifdef DIP_DIAGNOSTICS
    interpolate_calls++;
    flips += n_flips;
    avg_flips = flips / (float) interpolate_calls;
    
    // Welford's algorithm for calculating rolling mean and standard deviation
    double q = nn->get_quality();
    double delta = q - quality_avg;
    quality_avg += delta / interpolate_calls;
    M2_quality += delta * (q - quality_avg);
    
    if (interpolate_calls >= 2) {
      quality_stdev = M2_quality / (interpolate_calls - 1);
    }
    q_count++;
#endif
    

    return interp_vals;
}


int Cool::read_files(std::string cool_file, std::string tri_file, std::string neighbor_file) {
    /**
     * Reads files generated by CHIPS.
     *
     * cool_file        Path to file containing points of the grid (called *.points by default)
     * tri_file         Path to file containing triangulation (called *.tris by default)
     * neighbor_file   Path to file containing triangulation neighborhood relations (called *.neighbors by default)
     *
     * Returns 0 on success, otherwise 1, 2 or 3 if any of the files could not be read.
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
        std::cout << "Reading " << cool_file << std::endl;
    }

#ifdef DIP_POINTS_HEADER_SKIP
    std::getline(file, line);
#endif

    int n = 0;
    for (int i = 0; i < DIP_NMAX; i++) {
        std::getline(file, line);
        std::stringstream linestream(line);
        for (int j = 0; j < D; j++) {     // D coordinates
            std::getline(linestream, value, ',');
            points[i].coords[j] = std::stod(value);
        }
        for (int j = 0; j < DIP_VARNR; j++) {   // DIP_VARNR values
            std::getline(linestream, value, ',');
            points[i].value[j] = std::stod(value);
        }

        // Update smallest/largest known values
        // TODO Possibly replace with flags?
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
    N_LIM = n;

    file.close();


    /* Read Simplices */
    file.open(tri_file);
    if (!file.is_open()) {
        std::cerr << "Error reading " << tri_file << std::endl;
        return 2;
    } else {
        std::cout << "Reading " << tri_file << std::endl;
    }


    int s = 0;     // Counter to determine number of simplices
    int s_err = 0; // Number of errors in simplex loading process

    for (int i = 0; i < S_MAX; i++) {
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
        simplices[i].index = i;

        s++;
        int errbool = 0;
        if (skip != 1) {
            errbool += simplices[i].calculate_normals();
            errbool += simplices[i].validate_normals();
        }
        if (errbool) {
            s_err++;
        }

        // TODO Reinstate validate simplex function
//        simplices[i].validate_simplex();
        if (file.peek() == EOF) {
            break;
        }
    }
    S_LIM = s;


    file.close();

    std::cout << "DIP: " << s_err << " simplices out of " << S_LIM << " total had some kind of error. Does this match your expectations?" << std::endl;

    // Construct matrices for each simplex
    // https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Barycentric_coordinates_on_tetrahedra
    for (int s = 0; s < S_LIM; s++) {   // for each simplex
        simplices[s].construct_T_inv();
    }


    /* Read neighborhood relations */
    file.open(neighbor_file);
    if (!file.is_open()) {
        std::cerr << "Error reading " << neighbor_file << std::endl;
        return 3;
    } else {
        std::cout << "Reading " << neighbor_file << std::endl;
    }

    for (int i = 0; i < S_LIM; i++) {
        std::getline(file, line);
        std::stringstream linestream(line);
        for (int j = 0; j < D+1; j++) {     // D+1 points per simplex
            std::getline(linestream, value, ',');
            simplices[i].neighbor_indices[j] = std::stoi(value);

            if (std::stoi(value) != -1) {
                simplices[i].neighbor_pointers[j] = &(simplices[std::stoi(value)]);
            }
        }
    }

    file.close();

    return 0;
}


void Cool::set_clamp_values(double * cmins, double * cmaxs) {
    /**
     * Set values of the clamps to use in interpolate().
     *
     * double * cmins        Pointer to array of doubles. Min values for clamping.
     * double * cmaxs        Pointer to array of doubles. Max values for clamping.
     */
    for (int i = 0; i < DIP_DIMS; i++) {
        CLAMP_MAX[i] = cmaxs[i];
        CLAMP_MIN[i] = cmins[i];
    }
}