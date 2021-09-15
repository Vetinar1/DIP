//
// Created by vetinari on 01.12.20.
//

#ifndef DIP_COOLCONST_H
#define DIP_COOLCONST_H

#define DIP_NMAX 50000     // maximum number of points
//#define D 3             // dimensions
#define S_MAX 30000000    // maximum number of simplices


// The epsilon to use in DIP calculations
#define DIP_EPSILON 1e-6
#define AUTOMATIC_LOADING 1         // Whether to automatically load in new z slices in CoolManager. Warning: Assumes each call to interpolate has less or equal z than before
#define MAX_FLIPS 10                 // Max number of flips in simplex flipping algorithm before extrapolating

// The default value returned by interpolate if the ball tree fails
#define DIP_INTERP_DEFAULT (-40)

// Flag for asserts in ball tree generation code. Recommended to be kept on the first time a tree is built on a
// data set, can be switched off otherwise if desired
#define DIP_BALLTREE_CHECKS
// Whether to keep track of the number of interpolate calls, the total number of flips, and the average number of flips
#define DIP_DIAGNOSTICS
// Whether to skip the first line of the points file when reading files
#define DIP_POINTS_HEADER_SKIP

// Whether CoolManager objects should automatically load new z slices as required
#define DIP_CM_AUTOLOAD
// Whether CoolManager objects should quit the program upon encountering an error
//#define DIP_CM_ABORT_ON_ERROR

// Whether simplex errors in the file reading process are output to stderr or not
//#define DIP_SUPPRESS_SIMPLEX_ERRORS

#define DIP_DIMS 2
#define D DIP_DIMS
// The number of variables to be read/interpolated by DIP
// In order: Ctot, Htot, ne, H2, HI, HII, HeI, HeII, HeIII
#define DIP_VARNR 1

// #define PSI_SHOW_ERRORS
// #define PSI_SHOW_DIAGNOSTICS

#endif //DIP_COOLCONST_H
