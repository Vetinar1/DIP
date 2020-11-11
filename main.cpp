//
// Created by Stefan Lueders on 15.09.20.
//

#include <iostream>
#include <chrono>
#include "cool.cpp"
#include <climits>
#include <bitset>

int main() {
    const int N = 10;
    double grid[N][N][N];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                grid[i][j][k] = i*j*k;
            }
        }
    }

    double minmax[N][2] = {
            {0, 10},
            {0, 10},
            {0, 10}
    };

    int dims[3] = {10, 10, 10};

    MultilinearInterpolator<3> MLI(&grid[0][0][0], &minmax[0][0], &dims[0]);
    double point[3] = {2.33, 4.7, 3.1};
    double result = MLI.interpolate(&point[0]);
    std::cout << result << std::endl;

    return 0;
    std::cout << sizeof(Cool<221, 2, 436>) << std::endl;
    std::cout << "Initializing cool object... ";
//    Cool<35731, 4, 1002570> cool;
//    Cool<5994, 3, 38602> cool;
    Cool<1050, 2, 2082> cool;

    std::cout << "Done" << std::endl << "Reading files... ";
//    cool.read_files("../data4d/data.csv", "../data4d/dtri.csv", "../data4d/dneighbours.csv");
//    cool.read_files("../data3d/data.csv", "../data3d/dtri.csv", "../data3d/dneighbours.csv");
    cool.read_files("../data2d/data.csv", "../data2d/dtri.csv", "../data2d/dneighbours.csv");

    std::cout << "Done" << std::endl << "Constructing ball tree... " << std::flush;
    cool.construct_simplex_btree();
    std::cout << "Done" << std::endl << "Saving ball tree... ";
    cool.save_sbtree("tree");
    std::cout << "Done" << std::endl << "Beginning interpolation... " << std::endl;

    std::ofstream outfile;
    outfile.open("interp");

    double coord[2];
//    coord[2] = 0;
//    coord[3] = 20;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (int i = 0; i < 100; i++) {
        for (int j = 0; j < 100; j++) {
            std::cout << i << " " << j << std::endl;
            coord[0] = 2 + i * (8-2)/100.;
            coord[1] = -4 + j * 8 / 100.;
            double interp = cool.interpolate_sbtree(coord);

            outfile << coord[0] << " " << coord[1] << " " << interp << std::endl;
            std::cout << std::endl;
        }
    }
    std::cout << "Done" << std::endl;
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time to complete = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    std::cout << "Avg flips: " << cool.avg_sbtree_flips << std::endl;
    return 0;
}