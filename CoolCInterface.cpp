//
// Created by vetinari on 09.12.20.
//

#ifndef MASTER_PROJECT_C_PART_COOLCINTERFACE_CPP
#define MASTER_PROJECT_C_PART_COOLCINTERFACE_CPP

//#include "CoolCool.h"
#include "CoolCool.h"
#include <stdio.h>
#include <string>

// TODO Find a good solution for these

// Wrapper to be understood by C
// 1. Create or destroy Cool Objects
// 2. Handling of multiple objects?
// 3. Interface for member functions

//Cool<N, D, S> * C;
Cool * C;

#ifdef __cplusplus
extern "C" {
#endif

    void Cool_new() {
//        C = new Cool<N, D, S>;
        C = new Cool;
    }

    int Cool_read_files(char * cool_file, char * tri_file, char * neighbour_file) {
        return C->read_files(std::string(cool_file), std::string(tri_file), std::string(neighbour_file));
    }

    void Cool_reset() {
        C->reset();
    }

    void Cool_construct_tree() {
        C->construct_btree();
    }

    void Cool_save_tree(char* fname) {
        C->save_btree(std::string(fname));
    }

    double Cool_interpolate(double * coords) {
        return C->interpolate(coords);
    }


    /* CoolManager<int N_LIM, int D, int S_LIM> */
    // TODO Find way to pass filenames
//    CoolManager<N_LIM, D, S_LIM> * CoolManager_new(double init_z_low, double init_z_high, ) {
//        return new CoolManager<N_LIM, D, S_LIM>;
//    }
//
//    double CoolManager_interpolate(CoolManager<N_LIM, D, S_LIM> * CM, double * args, double z) {
//        return CM->interpolate(args, z);
//    }
//
//    void CoolManager_push_slice(CoolManager<N_LIM, D, S_LIM> * CM, char * fname) {
//        CM->push_slice(fname);
//    }

#ifdef __cplusplus
};
#endif


#endif //MASTER_PROJECT_C_PART_COOLCINTERFACE_CPP
