//
// Created by vetinari on 01.12.20.
//

#ifndef MASTER_PROJECT_C_PART_COOLCONST_H
#define MASTER_PROJECT_C_PART_COOLCONST_H

#define N_MAX 61000
#define D 4
#define S_MAX 1822576

//#define N_LIM 10000
//#define S_LIM 10000

#define EPSILON 1e-6
#define AUTOMATIC_LOADING 1         // Whether to automatically load in new z slices in CoolManager. Warning: Assumes each call to interpolate has less or equal z than before
#define MAX_FLIPS 5                 // Max number of flips in simplex flipping algorithm before extrapolating TODO log somehow?


// Flag for asserts in ball tree generation code. Recommended to be kept on the first time a tree is built on a
// data set, can be switched off otherwise if desired
#define DIP_BALLTREE_CHECKS
// Whether to keep track of the number of interpolate calls, the total number of flips, and the average number of flips
#define DIP_DIAGNOSTICS
// Whether to skip the first line of the points file when reading files
#define DIP_POINTS_HEADER_SKIP

#endif //MASTER_PROJECT_C_PART_COOLCONST_H
