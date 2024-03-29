//
// Created by vetinari on 26.11.20.
//

#ifndef DIP_COOLMANAGER_H
#define DIP_COOLMANAGER_H

#include "Cool.h"
#include <map>
#include <iostream>
#include <assert.h>
#include "Const.h"
#include "PSI.h"

class CoolManager {
    /**
     * Class that manages two Cool objects ("slices"), which are used for interpolation in all dimensions except z. z is
     * interpolated linearly between those.
     */
private:
#ifdef DIP_CM_USE_PSI
    PSI * low;
    PSI * high;
#else
    Cool * low;
    Cool * high;
#endif

    double z_low, z_high, z_diff, z_highest;
    double CLAMP_MIN[DIP_DIMS], CLAMP_MAX[DIP_DIMS];

    int highest_z_warn_flag = 0;
    int z_high_warning_val = 9999999;

    std::map<double, std::string> filenames;    // Available files for each z. (endings: .points, .tris, .neighbors)

    void autoload(double z);

    void apply_clamps();
public:
    CoolManager(double init_z_low, double init_z_high, std::string fname) {
        std::ifstream mapfile;
        mapfile.open(fname);

        if (!mapfile.is_open()) {
            std::cerr << "Error reading mapfile" << std::endl;
            abort();
        }

        double key;
        std::string value;
        while (mapfile >> key >> value) {
            filenames[key] = value;
//            std::cout << key << " " << value << " " << filenames[key] << std::endl;
        }
        mapfile.close();

        assert(init_z_high > init_z_low);
        z_low = init_z_low;
        z_high = init_z_high;
        z_diff = init_z_high - init_z_low;
        z_highest = z_high;

#ifdef DIP_CM_USE_PSI
        low  = new PSI;
        high = new PSI;
    
        int status = low->read_files(filenames[z_low] + ".points");
        if (status > 0) {
            std::cerr << "Error reading cooling data in CoolManager init, files " << filenames[z_low] << std::endl;
            abort();
        }
        status = high->read_files(filenames[z_high] + ".points");
        if (status > 0) {
            std::cerr << "Error reading cooling data in CoolManager init, files " << filenames[z_high] << std::endl;
            abort();
        }
#else
        low  = new Cool;
        high = new Cool;

        int status = low->read_files(
            filenames[z_low] + ".points",
            filenames[z_low] + ".tris",
            filenames[z_low] + ".neighbors"
        );
        if (status > 0) {
            std::cerr << "Error reading cooling data in CoolManager init, files " << filenames[z_low] << std::endl;
            abort();
        }
        status = high->read_files(
            filenames[z_high] + ".points",
            filenames[z_high] + ".tris",
            filenames[z_high] + ".neighbors"
        );
        if (status > 0) {
            std::cerr << "Error reading cooling data in CoolManager init, files " << filenames[z_high] << std::endl;
            abort();
        }

        low->construct_btree();
        high->construct_btree();
#endif
    }
    ~CoolManager() {
        delete low;
        delete high;
    }
    double * interpolate(double*, double);
    void push_slice(std::string);
    void set_clamp_values(double * mins, double * maxs);
};


#endif // DIP_COOLMANAGER_H