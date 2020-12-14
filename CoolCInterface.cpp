//
// Created by vetinari on 09.12.20.
//

#ifndef MASTER_PROJECT_C_PART_COOLCINTERFACE_CPP
#define MASTER_PROJECT_C_PART_COOLCINTERFACE_CPP

#include "CoolCool.h"
#include "CoolManager.h"

// TODO Find a good solution for these
#define N 10000
#define D 5
#define S 10000

#define N_MAX 10000
#define S_MAX 10000

// Wrapper to be understood by C
// 1. Create or destroy Cool Objects
// 2. Handling of multiple objects?
// 3. Interface for member functions

extern "C" {

    /* Cool<int N, int D, int S> */
    // TODO Can I return pointers like this? Probably not
    // Work with indices and a list?
    Cool<N, D, S> * Cool_new() {
        return new Cool<N, D, S>;
    }

    int Cool_read_files(Cool<N, D, S> * C, char * cool_file, char * tri_file, char * neighbour_file) {
        return C->read_files(std::string(cool_file), std::string(tri_file), std::string(neighbour_file));
    }

    void Cool_reset(Cool<N, D, S> * C) {
        C->reset();
    }

    void Cool_construct_tree(Cool<N, D, S> * C) {
        C->construct_btree();
    }

    void Cool_save_tree(Cool<N, D, S> * C, char* fname) {
        C->save_btree(std::string(fname));
    }

    double Cool_interpolate(Cool<N, D, S> * C, double * coords) {
        return C->interpolate(coords);
    }


    /* CoolManager<int N_MAX, int D, int S_MAX> */
    // TODO Find way to pass filenames
//    CoolManager<N_MAX, D, S_MAX> * CoolManager_new(double init_z_low, double init_z_high, ) {
//        return new CoolManager<N_MAX, D, S_MAX>;
//    }

    double CoolManager_interpolate(CoolManager<N_MAX, D, S_MAX> * CM, double * args, double z) {
        return CM->interpolate(args, z);
    }

    void CoolManager_push_slice(CoolManager<N_MAX, D, S_MAX> * CM, char * fname) {
        CM->push_slice(fname);
    }

};


// Undef to not interfere with implementation. TODO: Correct? Rename?
#undef N
#undef D
#undef S

#undef N_MAX
#undef S_MAX

#endif //MASTER_PROJECT_C_PART_COOLCINTERFACE_CPP
