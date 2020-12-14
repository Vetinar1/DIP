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
#include <unordered_map>
#include <algorithm>
#include <bitset>
#include "CoolPoint.h"
#include "CoolSimplex.h"
#include "CoolConst.h"


template<int N, int D, int S>
class Cool {
    /**
     * int N        Size of points array. Recommended to be number of lines in .points file. More will also work.
     * int D        Dimensionality of parameter space.
     * int S        Size of simplices array. Recommended to be number of lines .tris/.neighbors files. More also works.
     *
     * This is the main object for reading and interpolating cooling data.
     * Start by reading in the .points, .tris, and .neighbors files generated by CHIPS using read_files(). Then
     * construct the ball tree using construct_btree(). The object is then ready for interpolation using the
     * interpolate() method.
     * In general, coordinates (e.g. to interpolate) are passed as pointers to arrays. These arrays need to have correct
     * length (matching dimension D) and contain the coordinates in the same order as the CHIPS data files.
     * If you want to read in new data, but don't want to allocate a new Cool object, you can use reset().
     *
     * Make sure to allocate Cool objets on the heap, otherwise the size of the points and simplices attributes will
     * lead to a stack overflow.
     */
private:

    Point<D> points[N];
    Simplex<D> simplices[S];
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

    // Minimum and maximum values in each dimension TODO private?
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



#endif //MASTER_PROJECT_C_PART_COOLCOOL_H