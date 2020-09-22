//
// Created by Stefan Lueders on 15.09.20.
//

#include <iostream>
#include <chrono>
#include "cool.cpp"

int main() {
    Cool<221, 2, 436> cool;

    cool.read_files("../data.csv", "../dtri.csv", "../dneighbours.csv");
    cool.construct_btree();
    cool.save_tree("tree");

    std::ofstream outfile;
    outfile.open("interp");

    double coord[2];
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (int i = 0; i < 100; i++) {
        for (int j = 0; j < 100; j++) {
            coord[0] = 2 + i * (8-2)/100.;
            coord[1] = -4 + j * 8 / 100.;
            double interp = cool.interpolate(coord);

            outfile << coord[0] << " " << coord[1] << " " << interp << std::endl;
        }
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time to complete = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    std::cout << "Avg flips: " << cool.avg_flips << std::endl;
    return 0;
}