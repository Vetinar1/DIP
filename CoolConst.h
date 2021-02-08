//
// Created by vetinari on 01.12.20.
//

#ifndef MASTER_PROJECT_C_PART_COOLCONST_H
#define MASTER_PROJECT_C_PART_COOLCONST_H

#define N 57076
#define D 4
#define S 1706632

//#define N_LIM 10000
//#define S_LIM 10000

#define EPSILON 1e-6
#define AUTOMATIC_LOADING 1         // Whether to automatically load in new z slices in CoolManager. Warning: Assumes each call to interpolate has less or equal z than before
#define MAX_FLIPS 100               // Max number of flips in simplex flipping algorithm before throwing an error
#define COOL_CLAMP_COORDS 1         // Whether to automatically input coordinates to interpolate call. Warning: If enabled may modify inputs in place!

#endif //MASTER_PROJECT_C_PART_COOLCONST_H
