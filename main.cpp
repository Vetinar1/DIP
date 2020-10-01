//
// Created by Stefan Lueders on 15.09.20.
//

#include <iostream>
#include <chrono>
#include "cool.cpp"

int main() {
//    std::cout << sizeof(double[35090][4]) / 1e6 << std::endl;
//    std::cout << sizeof(Simplex<4>[1002570]) / 1e6 << std::endl;
    std::cout << "Initializing cool object... ";
    Cool<35731, 4, 1002570> cool;

    std::cout << "Done" << std::endl << "Reading files... ";
    cool.read_files("../data4d/data.csv", "../data4d/dtri.csv", "../data4d/dneighbours.csv");
    std::cout << "Done" << std::endl << "Constructing ball tree... ";
    cool.construct_simplex_btree();
//    std::cout << "Done" << std::endl << "Saving ball tree... ";
//    cool.save_sbtree("tree");
    std::cout << "Done" << std::endl << "Beginning interpolation... " << std::endl;

    std::ofstream outfile;
    outfile.open("interp");

    double coord[4];
    coord[2] = 20;
    coord[3] = 20;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (int i = 0; i < 100; i++) {
        for (int j = 0; j < 100; j++) {
            std::cout << i << " " << j << std::endl;
            coord[0] = 2 + i * (8-2)/100.;
            coord[1] = -4 + j * 8 / 100.;
            double interp = cool.interpolate_sbtree(coord);

            outfile << coord[0] << " " << coord[1] << " " << interp << std::endl;
        }
    }
    std::cout << "Done" << std::endl;
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time to complete = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    std::cout << "Avg flips: " << cool.avg_sbtree_flips << std::endl;
    return 0;
}