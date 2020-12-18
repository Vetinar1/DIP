//
// Created by vetinari on 26.11.20.
//

#ifndef MASTER_PROJECT_C_PART_COOLMANAGER_H
#define MASTER_PROJECT_C_PART_COOLMANAGER_H

#include "CoolCool.h"
#include <map>
#include <iostream>
#include <assert.h>
#include "CoolConst.h"

class CoolManager {
    /**
     * Class that manages Cool objects in order to do linear interpolation between z slices.
     * TODO better description
     */
private:
    Cool * low;
    Cool * high;

    double z_low, z_high, z_diff;

    std::map<double, std::string> filenames;    // Available files for each z. (endings: .points, .tris, .neighbors)

    void autoload(double z);
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

        low  = new Cool;
        high = new Cool;

        low->read_files(
                filenames[z_low] + ".points",
                filenames[z_low] + ".tris",
                filenames[z_low] + ".neighbors"
        );
        high->read_files(
                filenames[z_high] + ".points",
                filenames[z_high] + ".tris",
                filenames[z_high] + ".neighbors"
        );

        low->construct_btree();
        high->construct_btree();
    }
    ~CoolManager() {
        delete low;
        delete high;
    }
    double interpolate(double*, double);
    void save_trees(std::string fname_low, std::string fname_high);
    void push_slice(std::string);
};


#endif // MASTER_PROJECT_C_PART_COOLMANAGER_H