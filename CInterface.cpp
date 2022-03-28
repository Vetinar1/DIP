/**
 * This file makes methods of Cool, CoolManager and PSI objects available in C.
 * Please check the respective method implementation for documentation.
 */

#ifndef DIP_CInterface_CPP
#define DIP_CInterface_CPP

#include "Cool.h"
#include "CoolManager.h"
#include <stdio.h>
#include <string>
#include "CInterface.h"
#include "PSI.h"

#define COOL_OBJ_COUNT 10
#define COOLMANAGER_OBJ_COUNT 161
#define PSI_OBJ_COUNT 10

Cool * Cs[COOL_OBJ_COUNT];
CoolManager * CMs[COOLMANAGER_OBJ_COUNT];
PSI * Ps[PSI_OBJ_COUNT];

int cool_counter = 0;
int coolmanager_counter = 0;
int psi_counter = 0;

#ifdef __cplusplus
extern "C" {
#endif

    int Cool_new(int c_indx) {
        if (cool_counter >= COOL_OBJ_COUNT) {
            std::cerr << "Tried to create too many Cool objects" << std::endl;
            return -1;
        }
        Cs[c_indx] = new Cool;
        cool_counter++;
        return cool_counter - 1;
    }

    int Cool_read_files(int c_indx, char * cool_file, char * tri_file, char * neighbor_file) {
        return Cs[c_indx]->read_files(std::string(cool_file), std::string(tri_file), std::string(neighbor_file));
    }

    void Cool_reset(int c_indx) {
        Cs[c_indx]->reset();
    }

    void Cool_construct_tree(int c_indx) {
        Cs[c_indx]->construct_btree();
    }

    void Cool_save_tree(int c_indx, char* fname) {
        Cs[c_indx]->save_btree(std::string(fname));
    }

    double * Cool_interpolate(int c_indx, double * coords) {
        return Cs[c_indx]->interpolate(coords);
    }

    void Cool_set_clamps(int c_indx, double * mins, double * maxs) {
        Cs[c_indx]->set_clamp_values(mins, maxs);
    }


    int CoolManager_new(double init_z_low, double init_z_high, char * mapfile, int index) {
        if (coolmanager_counter >= COOLMANAGER_OBJ_COUNT) {
            std::cerr << "Tried to create too many CoolManager objects" << std::endl;
            return -1;
        }
        CMs[index] = new CoolManager(init_z_low, init_z_high, std::string(mapfile));
        coolmanager_counter++;
        return index;
    }

    double * CoolManager_interpolate(int cm_indx, double * args, double z) {
        return CMs[cm_indx]->interpolate(args, z);
    }

    void CoolManager_set_clamps(int cm_idx, double * mins, double * maxs) {
        CMs[cm_idx]->set_clamp_values(mins, maxs);
    }

    
    int PSI_new(int p_indx) {
        if (psi_counter >= PSI_OBJ_COUNT) {
            std::cerr << "Tried to create too many PSI objects" << std::endl;
            return -1;
        }
        Ps[p_indx] = new PSI;
        psi_counter++;
        return psi_counter - 1;
    }
    
    int PSI_read_files(int p_indx, char * cool_file) {
        return Ps[p_indx]->read_files(std::string(cool_file));
    }
    
    void PSI_reset(int p_indx) {
        Ps[p_indx]->reset();
    }
    
    double * PSI_interpolate(int p_indx, double * coords, int smart_fallback) {
        return Ps[p_indx]->interpolate(coords, smart_fallback);
    }
    
    void PSI_set_clamps(int p_indx, double * mins, double * maxs) {
        Ps[p_indx]->set_clamp_values(mins, maxs);
    }

#ifdef __cplusplus
};
#endif


#endif //DIP_CInterface_CPP
