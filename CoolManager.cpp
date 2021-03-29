//
// Created by vetinari on 14.12.20.
//

#include "CoolCool.h"
#include <map>
#include <iostream>
#include <assert.h>
#include "CoolConst.h"
#include "CoolManager.h"


double CoolManager::interpolate(double * args, double z) {
    if (z < z_low || z > z_high) {
        if (z > z_high && z_high == z_highest) {
            if (highest_z_warn_flag != 1) {
                std::cerr << "DIP Warning! Input z greater than largest z (" << z_highest << ")." << std::endl;
                std::cerr << "Assuming that data at z = " << z_high << " generalizes to " << z << std::endl;
                std::cerr << "This warning will only be output once" << std::endl;
                highest_z_warn_flag = 1;
            }
            return high->interpolate(args);
        }
#ifdef DIP_CM_AUTOLOAD
        else {
            autoload(z);
        }
#endif
    }
    // Interpolate inside slices
    double lambda_low = low->interpolate(args);
//    std::cout << "Low: " << lambda_low << std::endl;
    double lambda_high = high->interpolate(args);
//    std::cout << "High: " << lambda_high << std::endl;
    // Simple linear interpolation in between slices
    return (lambda_low * (z_high - z) + lambda_high * (z - z_low)) / z_diff;
}


void CoolManager::autoload(double z) {
    /**
     * If the current z is outside of the bounds of the current two slices, load another slice so that z is covered.
     * Note: Only updates one slice. This means it implicitly assumes we only ever move to adjacent slices, forwards
     * in time! This should be a reasonable assumption for a cosmological simulation.
     */
    if (z <= z_high && z >= z_low) {
        return;
    }

    if (z > z_high) {
        if (z_high < z_high_warning_val) {
            z_high_warning_val = z_high;
            std::cerr << "DIP WARNING in CoolManager object " << this << " in function void CoolManager::autoload" << std::endl;
            std::cerr << "Moving backwards in time?" << std::endl;
            std::cerr << "Current z interval: [" << z_low << ", " << z_high << "]" << std::endl;
            std::cerr << "Tried to move to z = " << z << std::endl;
            std::cerr << "This warning will only be output once per z" << std::endl;
#ifdef DIP_CM_ABORT_ON_ERROR
            std::cerr << "Aborting..." << std::endl;
            abort();
#endif
        }
        return;
    }

    std::string fname;
    // std::map is sorted by keys
    for (auto it=filenames.begin(); it!=filenames.end(); ++it) {
        if (it->first >= z_low) {
            it--;
            fname = it->second;

            z_high = z_low;
            z_low = it->first;
            z_diff = z_high - z_low;
            break;
        }
    }

    push_slice(fname);
    apply_clamps();
}


void CoolManager::push_slice(std::string filename) {
    // TODO Manual usage of this function needs to update z_low, z_high and z_diff
    Cool * temp = high;
    high = low;
    low = temp;
    low->reset();
    int status = low->read_files(
            filename + ".points",
            filename + ".tris",
            filename + ".neighbors"
    );
    if (status > 0) {
        std::cerr << "Error reading cooling data in CoolManager::push_slice, files " << filename << std::endl;
        abort();
    }
    low->construct_btree();
}


void CoolManager::save_trees(std::string fname_low, std::string fname_high) {
    low->save_btree(fname_low);
    high->save_btree(fname_high);
}


void CoolManager::set_clamp_values(double * cmins, double * cmaxs) {
    /**
     * Set values of the clamps to use in interpolate().
     *
     * double * cmins        Pointer to array of doubles. Min values for clamping.
     * double * cmaxs        Pointer to array of doubles. Max values for clamping.
     */
    for (int i = 0; i < D; i++) {
        CLAMP_MAX[i] = cmaxs[i];
        CLAMP_MIN[i] = cmins[i];
    }

    apply_clamps();
}


void CoolManager::apply_clamps() {
    low->set_clamp_values(&CLAMP_MIN[0], &CLAMP_MAX[0]);
    high->set_clamp_values(&CLAMP_MIN[0], &CLAMP_MAX[0]);
}