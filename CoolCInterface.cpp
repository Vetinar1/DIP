//
// Created by vetinari on 09.12.20.
//

#ifndef MASTER_PROJECT_C_PART_COOLCINTERFACE_CPP
#define MASTER_PROJECT_C_PART_COOLCINTERFACE_CPP

#include "CoolCool.h"
#include "CoolManager.h"
#include <stdio.h>
#include <string>
//#include "CoolCInterface.h"

#define COOL_OBJ_COUNT 10
#define COOLMANAGER_OBJ_COUNT 161

// TODO Find a good solution for these

// Wrapper to be understood by C
// 1. Create or destroy Cool Objects
// 2. Handling of multiple objects?
// 3. Interface for member functions

//Cool<N, D, S> * C;
Cool * C[COOL_OBJ_COUNT];
CoolManager * CM[COOLMANAGER_OBJ_COUNT];

int cool_counter = 0;
int coolmanager_counter = 0;

#ifdef __cplusplus
extern "C" {
#endif

    int Cool_new() {
//        C = new Cool<N, D, S>;
        if (cool_counter >= COOL_OBJ_COUNT) {
            std::cerr << "Tried to create too many Cool objects" << std::endl;
            return -1;
        }
        C[cool_counter] = new Cool;
        cool_counter++;
        return cool_counter - 1;    // Less pretty but also less ambiguous
    }

    int Cool_read_files(int c_indx, char * cool_file, char * tri_file, char * neighbour_file) {
        return C[c_indx]->read_files(std::string(cool_file), std::string(tri_file), std::string(neighbour_file));
    }

    void Cool_reset(int c_indx) {
        C[c_indx]->reset();
    }

    void Cool_construct_tree(int c_indx) {
        C[c_indx]->construct_btree();
    }

    void Cool_save_tree(int c_indx, char* fname) {
        C[c_indx]->save_btree(std::string(fname));
    }

    double Cool_interpolate(int c_indx, double * coords) {
        return C[c_indx]->interpolate(coords);
    }


    /* CoolManager<int N_LIM, int D, int S_LIM> */
//    int CoolManager_new(double init_z_low, double init_z_high, char * mapfile) {
//        if (coolmanager_counter >= COOLMANAGER_OBJ_COUNT) {
//            std::cerr << "Tried to create too many CoolManager objects" << std::endl;
//            return -1;
//        }
//        CM[coolmanager_counter] = new CoolManager(init_z_low, init_z_high, std::string(mapfile));
//        coolmanager_counter++;
//        return coolmanager_counter - 1;
//    }

    int CoolManager_new(double init_z_low, double init_z_high, char * mapfile, int index) {
        /**
         * Dont mix this one with the other one
         */
        if (coolmanager_counter >= COOLMANAGER_OBJ_COUNT) {
            std::cerr << "Tried to create too many CoolManager objects" << std::endl;
            return -1;
        }
        CM[index] = new CoolManager(init_z_low, init_z_high, std::string(mapfile));
        coolmanager_counter++;
        return index;
    }

    double CoolManager_interpolate(int cm_indx, double * args, double z) {
        return CM[cm_indx]->interpolate(args, z);
    }

    void CoolManager_save_trees(int cm_indx, char * fname_low, char * fname_high) {
        CM[cm_indx]->save_trees(std::string(fname_low), std::string(fname_high));
    }

    void CoolManager_set_clamps(int cm_idx, double * mins, double * maxs) {
        CM[cm_idx]->set_clamp_values(mins, maxs);
    }

#ifdef __cplusplus
};
#endif


#endif //MASTER_PROJECT_C_PART_COOLCINTERFACE_CPP
