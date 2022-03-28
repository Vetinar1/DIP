//
// Created by vetinari on 14.12.20.
//

#include "Cool.h"
#include <map>
#include <iostream>
#include <assert.h>
#include "Const.h"
#include "CoolManager.h"


double * CoolManager::interpolate(double * args, double z) {
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
#if defined(DIP_CM_USE_PSI) && defined(DIP_CM_USE_PSI_FALLBACK)
    double * low_vals = low->interpolate(args, 1);
    double * high_vals = high->interpolate(args, 1);
#elif defined(DIP_CM_USE_PSI)
    double * low_vals = low->interpolate(args, 0);
    double * high_vals = high->interpolate(args, 0);
#else
    double * low_vals = low->interpolate(args);
    double * high_vals = high->interpolate(args);
#endif
    
    if (low_vals == nullptr && high_vals == nullptr) {
        return nullptr;
    } else if (low_vals == nullptr) {
        return high_vals;
    } else if (high_vals == nullptr) {
        return low_vals;
    }
    
    double * interp_vals = new double[DIP_VARNR];
    // Simple linear interpolation in between slices
    for (int i = 0; i < DIP_VARNR; i++) {
        interp_vals[i] = (low_vals[i] * (z_high - z) + high_vals[i] * (z - z_low)) / z_diff;
    }
    
    delete[] low_vals;
    delete[] high_vals;
    
//    return (lambda_low * (z_high - z) + lambda_high * (z - z_low)) / z_diff;
    return interp_vals;
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
    /**
     * Unload upper slice, push lower to upper slice, add new lower slice
     *
     * Manual usage of this function needs to update z_low, z_high and z_diff
     */
#ifdef DIP_CM_USE_PSI
    PSI * temp = high;
    high = low;
    low = temp;
    low->reset();
    int status = low->read_files(filename + ".points");
#else
    Cool * temp = high;
    high = low;
    low = temp;
    low->reset();
    int status = low->read_files(
        filename + ".points",
        filename + ".tris",
        filename + ".neighbors"
    );
    if (status == 0) {
        low->construct_btree();
    }
#endif
    if (status > 0) {
        std::cerr << "Error reading cooling data in CoolManager::push_slice, files " << filename << std::endl;
        abort();
    }
}


void CoolManager::set_clamp_values(double * cmins, double * cmaxs) {
    /**
     * Set values of the clamps to use in interpolate().
     *
     * double * cmins        Pointer to array of doubles. Min values for clamping.
     * double * cmaxs        Pointer to array of doubles. Max values for clamping.
     */
    for (int i = 0; i < DIP_DIMS; i++) {
        CLAMP_MAX[i] = cmaxs[i];
        CLAMP_MIN[i] = cmins[i];
    }

    apply_clamps();
}


void CoolManager::apply_clamps() {
    low->set_clamp_values(&CLAMP_MIN[0], &CLAMP_MAX[0]);
    high->set_clamp_values(&CLAMP_MIN[0], &CLAMP_MAX[0]);
}