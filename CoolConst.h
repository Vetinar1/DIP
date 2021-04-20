//
// Created by vetinari on 01.12.20.
//

#ifndef MASTER_PROJECT_C_PART_COOLCONST_H
#define MASTER_PROJECT_C_PART_COOLCONST_H

#define N_MAX 14000     // maximum number of points
#define D 4             // dimensions
#define S_MAX 410000    // maximum number of simplices


// The epsilon to use in DIP calculations
#define DIP_EPSILON 1e-6
#define AUTOMATIC_LOADING 1         // Whether to automatically load in new z slices in CoolManager. Warning: Assumes each call to interpolate has less or equal z than before
#define MAX_FLIPS 10                 // Max number of flips in simplex flipping algorithm before extrapolating


// Flag for asserts in ball tree generation code. Recommended to be kept on the first time a tree is built on a
// data set, can be switched off otherwise if desired
#define DIP_BALLTREE_CHECKS
// Whether to keep track of the number of interpolate calls, the total number of flips, and the average number of flips
#define DIP_DIAGNOSTICS
// Whether to skip the first line of the points file when reading files
#define DIP_POINTS_HEADER_SKIP

// Starting values for CoolManager
#define DIP_CM_INIT_Z_LOW 8.5
#define DIP_CM_INIT_Z_HIGH 9.
#define DIP_CM_MAPFILE "cooldata/mapfile"

// Whether CoolManager objects should automatically load new z slices as required
#define DIP_CM_AUTOLOAD
// Whether CoolManager objects should quit the program upon encountering an error
//#define DIP_CM_ABORT_ON_ERROR

// Whether simplex errors in the file reading process are output to stderr or not
//#define DIP_SUPPRESS_SIMPLEX_ERRORS

#endif //MASTER_PROJECT_C_PART_COOLCONST_H
