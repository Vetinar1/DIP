//
// Created by vetinari on 01.12.20.
//

#ifndef DIP_COOLCONST_H
#define DIP_COOLCONST_H

#define DIP_NMAX 50000     // maximum number of points
//#define D 3             // dimensions
#define S_MAX 41000000    // maximum number of simplices


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

// Whether CoolManager objects should automatically load new z slices as required
#define DIP_CM_AUTOLOAD
// Whether CoolManager objects should quit the program upon encountering an error
//#define DIP_CM_ABORT_ON_ERROR
// Whether CoolManager objects should use PSI objects instead of Cool objects
//#define DIP_CM_USE_PSI

// Whether simplex errors in the file reading process are output to stderr or not
//#define DIP_SUPPRESS_SIMPLEX_ERRORS

#define DIP_DIMS 3
#define D DIP_DIMS
// The number of variables to be read/interpolated by DIP
// In order: Htot, Ctot, ne, H2, HI, HII, HeI, HeII, HeIII
#define DIP_VARNR 9

// #define PSI_SHOW_ERRORS
// #define PSI_SHOW_DIAGNOSTICS



// Constants specifically for PSI
// Number of nearest neighbors to use in PSI algorithm
#define PSI_K 50
// If simplex construction fails, multiply PSI_K by this and try again
#define PSI_KFACTOR 2
// Maximum number of failed simplex constructions before using default value (DIP_INTERP_DEFAULT)
#define PSI_MAXREP 4

#endif //DIP_COOLCONST_H
