//
// Created by vetinari on 01.12.20.
//

#ifndef MASTER_PROJECT_C_PART_COOLMLI_H
#define MASTER_PROJECT_C_PART_COOLMLI_H

#include <math.h>
#include <bitset>

//template<int N>
class MultilinearInterpolator {
private:
    int dims[N];
    double * grid;
    double minmax[N][2];
    double dim_lens[N];
    double dim_offsets[N];

    double linear_interpolation(double, double, double, double, double);
public:
    MultilinearInterpolator(double * in_grid, const double * in_minmax, const int * in_dims) {
        grid = in_grid;
        for (int i = 0; i < N; i++) {
            minmax[i][0] = in_minmax[2*i];
            minmax[i][1] = in_minmax[2*i+1];
            dim_lens[i] = minmax[i][1] - minmax[i][0];
        }
        for (int i = 0; i < N; i++) {
            dims[i] = in_dims[i];
        }
        for (int i = 0; i < N; i++) {
            dim_offsets[i] = 1;
            if (i == N) {
                break;
            }
            for (int j = i+1; j < N; j++) {
                dim_offsets[i] *= dims[j];
            }
        }
    }
    double interpolate(const double *);
};


//template<int N>
double MultilinearInterpolator<N>::linear_interpolation(double x, double x0, double y0, double x1, double y1) {
    return (y0 * (x1 - x) + y1 * (x - x0)) / (x1 - x0);
}

//template<int N>
double MultilinearInterpolator<N>::interpolate(const double* point) {
    /**
     * Optimized multilinear interpolation function.
     * Undefined behaviour when trying to interpolate outside grid.
     *
     * point = array of length N, point to interpolate at
     */

    // 1. Figure out points of surrounding hypercube
    int indices[N];
    double diffs[N];

    for (int i = 0; i < N; i++) {
        // Automatically cast to int!
        indices[i] = (point[i] - minmax[i][0]) * dims[i] / dim_lens[i];
        double xi_0 = minmax[i][0] + dim_lens[i] * indices[i] / dims[i];
        double xi_1 = minmax[i][0] + dim_lens[i] * (indices[i] + 1) / dims[i];
        diffs[i] = (point[i] - xi_0) / (xi_1 - xi_0);
    }

    double * vals = new double[(int) pow(2, N)];
    for (int i = 0; i < pow(2, N); i++) {
        std::bitset<N> offsets(i);
        int index = 0;

        for (int j = 0; j < N; j++) {
            index += (offsets[j]) ? dim_offsets[j] * (indices[j] + 1) : dim_offsets[j] * indices[j];
        }

        vals[i] = grid[index];
    }

    double * new_vals;
    for (int i = N; i > 0; i--) {
        new_vals = new double[(int) pow(2, i-1)];
        for (int j = 0; j < (int) pow(2, i-1); j++) {
            // To understand why the offset in the second vals[] works, look at the bit representation of the indices
            // Those bits correspond to the offsets above
            new_vals[j] = vals[j] * (1 - diffs[i]) + vals[j + (int) pow(2, i-1)] * diffs[i];
        }

        delete[] vals;
        vals = new_vals;
    }

    double result = new_vals[0];
    delete[] new_vals;
    return result;
}

#endif //MASTER_PROJECT_C_PART_COOLMLI_H
