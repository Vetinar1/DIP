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
    int points[D+1];                    // D+1 points (indices)
    int neighbour_indices[D+1];         // One neighbour oppposite every point
    Simplex * neighbour_pointers[D+1];
    double centroid[D];
    double btree_radius_sq;
public:
    Simplex * lchild;
    Simplex * rchild;

    Simplex() {
        btree_radius_sq = 0;
    };

    void calculate_centroid(const double * coords) {
        /**
         * double * coords      Pointer to array to use for lookup of actual coordinates - i.e. Cool.points
         */
        for (int i = 0; i < D; i++) {   // coordinates
            centroid[i] = 0;
            for (int j = 0; j < D+1; j++) {     // points
                centroid[i] += *(coords + points[j]*(D+1) + i);
            }
            centroid[i] /= (D+1);
        }
    }
};

template<int N, int D, int S>
class Cool {
    /**
     * int N        Number of points
     * int D        Dimensionality of points
     * int S        Number of Simplices
     */
private:
    double points[N][D+1];      // N points, D dimensions, 1 Value
    Simplex<D> simplices[S];
    Simplex<D> * btree;         // Points to root of the ball tree

    Simplex<D> * construct_btree_recursive(Simplex<D> **, int, int);

public:
    Cool() {};

    int read_files(std::string, std::string, std::string);
    int construct_btree();
};

template<int N, int D, int S>
int Cool<N, D, S>::construct_btree() {
    /**
     * This function serves as a public "adapter" to the actual ball tree construction function,
     * construct_btree_recursive().
     *
     * Returns 0 on success.
     */
     // Construct array of pointers
     Simplex<D> * simps[S];
     for (int i = 0; i < S; i++) {
         simps[i] = &(simplices[i]);
     }

     printf("%i\n", S);
     btree = construct_btree_recursive(simps, S, 1);

     return 0;
}

template< int N, int D, int S>
Simplex<D> * Cool<N, D, S>::construct_btree_recursive(Simplex<D> ** simps, int n, int depth) {
    /**
     * Simplex<D> simps     Pointer to array of pointers to simplices to organize in tree
     * int n                Number of elements in array
     */
     printf("Recursion level: %i, n = %i\n", depth, n);

    // Input validation + what if theres only one element left?
    if (n < 0) {
        std::cerr << "Illegal argument in construct_btree: n = " << n << std::endl;
    } else if (n == 1) {
        simps[0]->lchild = NULL;
        simps[0]->rchild = NULL;
        simps[0]->btree_radius_sq = 0;

        printf("Recursing back from %i\n", depth);
        return *simps;
    }

    // More than one element
    // 1. Find dimension of greatest spread
    // Working with the centroids of the triangles
    double largest_spread = -1;
    int lspread_dim = -1;
    double avg = 0;             // Traditionally Ball trees work with the median, but that requires sorting

    for (int i = 0; i < D; i++) {   // every dimension
        double dim_min =    DBL_MAX;
        double dim_max = -1*DBL_MAX;
        double dim_spread = 0;
        for (int j = 0; j < n; j++) {   // every simplex
            double val = simps[j]->centroid[i];
            if (val < dim_min) {
                dim_min = val;
            } else if (val > dim_max) {
                dim_max = val;
            }
        }

        dim_spread = fabs(dim_max - dim_min);   // TODO: Should probably be some kind of assert instead of fabs...
        if (dim_spread > largest_spread) {
            largest_spread = dim_spread;
            lspread_dim    = i;
            avg            = (dim_max + dim_min) / 2;
        }
    }

    assert(largest_spread != -1 && lspread_dim != -1);

    // 2. Find pivot - simplex with centroid closest to avg
    // 3. Group points into sets to the left and right of average (L, R)
    int min_dist = DBL_MAX;
    int pivot_index;
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
            pivot_index = i;
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

    // 4. Recurse on L and R
    if (lcount > 0) {
        pivot_addr->lchild = construct_btree_recursive(L, lcount, depth+1);
    } else {
        pivot_addr->lchild = NULL;
    }
    if (rcount > 0) {
        pivot_addr->rchild = construct_btree_recursive(R, rcount, depth+1);
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
            dist += pow(simps[i]->centroid[j], 2);
        }

        if (dist > pivot_addr->btree_radius_sq) {
            pivot_addr->btree_radius_sq = dist;
        }
    }

    delete[] L;
    delete[] R;

    printf("Recursing back from %i\n", depth);

    return pivot_addr;
}

template<int N, int D, int S>
int Cool<N, D, S>::read_files(std::string cool_file, std::string tri_file, std::string neighbour_file) {
    /**
     * Reads files generated by the python program.
     *
     * cool_file        Path to file containing points of the grid
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
            this->points[i][j] = std::stod(value);
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
        for (int j = 0; j < D+1; j++) {     // D+1 points per simplex
            std::getline(linestream, value, ',');
            simplices[i].points[j] = std::stoi(value);
        }
        simplices[i].calculate_centroid(&(this->points[0][0]));
    }

    file.close();

    /* Read neighbourhood relations */
    file.open(neighbour_file);
    if (!file.is_open()) {
        std::cerr << "Error reading " << neighbour_file << std::endl;
        return 2;
    }


    for (int i = 0; i < S; i++) {
        std::getline(file, line);
        std::stringstream linestream(line);
        for (int j = 0; j < D+1; j++) {     // D+1 points per simplex
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