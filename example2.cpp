//
// Created by Stefan Lueders on 15.09.20.
//

#include <iostream>
#include <chrono>
#include "CoolCool.h"
#include "CoolManager.h"
#include <climits>
#include <bitset>
#include <math.h>
#include <string>
#include <map>
#include <iostream>

#include <iomanip>
#include <sstream>


int main() {
    double coord[3]; // T, nH, Z
    double z = 3.42;

    double clamp_mins[3] = {2.2, -8.8, -2.8};
    double clamp_maxs[3] = {8.8, 3.8, 0.8};

    std::cout << "Creating CoolManager object.." << std::endl;
    CoolManager * CM = new CoolManager(3.0, 3.5, "example_data/mapfile");

    std::cout << "Setting clamps..." << std::endl;
    CM->set_clamp_values(clamp_mins, clamp_maxs);

    std::cout << "Interpolating 10000 points using CoolManager" << std::endl;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            for (int k = 0; k < 10; k++) {
                coord[0] = clamp_mins[0] + i * (clamp_maxs[0] - clamp_mins[0]) / 10.;
                coord[1] = clamp_mins[1] + j * (clamp_maxs[1] - clamp_mins[1]) / 10.;
                coord[2] = clamp_mins[2] + k * (clamp_maxs[2] - clamp_mins[2]) / 10.;
                double interp = CM->interpolate(coord, z);
            }
        }
    }
    std::cout << "Done" << std::endl;
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time to complete = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;


    std::cout << "Creating Cool object... " << std::endl;
    Cool * cool = new Cool;

    cool->set_clamp_values(clamp_mins, clamp_maxs);

    std::cout << "Reading files... ";


    cool->read_files(
            "example_data/z3.0.points",
            "example_data/z3.0.tris",
            "example_data/z3.0.neighbors"
    );

    std::cout << "Done" << std::endl << "Constructing ball tree... " << std::flush;
    cool->construct_btree();

    std::cout << "Done" << std::endl << "Beginning interpolation... " << std::endl;
    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < 100; i++) {
        for (int j = 0; j < 100; j++) {
            coord[0] = clamp_mins[0] + i * (clamp_maxs[0]-clamp_mins[0])/100.;
            coord[1] = clamp_mins[1] + j * (clamp_maxs[1]-clamp_mins[1])/100.;
            double interp = cool->interpolate(coord);
        }
    }
    std::cout << "Done" << std::endl;
    end = std::chrono::steady_clock::now();
    std::cout << "Time to complete = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    std::cout << "Avg = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / pow(100, D) << "[ms]" << std::endl;
    std::cout << "Avg flips: " << cool->avg_flips << std::endl;
    return 0;
}
