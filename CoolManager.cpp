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
#if AUTOMATIC_LOADING==1
    autoload(z);
#endif
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
        std::cerr << "Moving backwards in time?" << std::endl;
        std::cerr << "Current z interval: [" << z_low << ", " << z_high << "]" << std::endl;
        std::cerr << "Tried to move to z = " << z << std::endl; // TODO This will lead to rpoblems if one calculation is ahead of the others
        assert(0);
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
}


void CoolManager::push_slice(std::string filename) {
    // TODO Manual usage of this function needs to update z_low, z_high and z_diff
    Cool * temp = high;
    high = low;
    low = temp;
    low->reset();
    low->read_files(
            filename + ".points",
            filename + ".tris",
            filename + ".neighbors"
    );
    low->construct_btree();
}


void CoolManager::save_trees(std::string fname_low, std::string fname_high) {
    low->save_btree(fname_low);
    high->save_btree(fname_high);
}