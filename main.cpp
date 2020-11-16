//
// Created by Stefan Lueders on 15.09.20.
//

#include <iostream>
#include <chrono>
#include "cool.cpp"
#include <climits>
#include <bitset>
#include <math.h>

int main() {
    const int X = 86;
    const int Y = 71;
    const int Z = 10;
    int dims[3] = {X, Y, Z};
    double grid[X][Y][Z];

    double minmax[3][2] = {
            {1, 9.5},
            {-8, 6},
            {-4, 0.5}
    };


    std::ifstream file;
    std::string line;
    std::string value;
    file.open("../repr1_schaye.csv");

    int x = 0;
    int y = 0;
    int z = 0;
    for (int i = 0; i < X*Y*Z; i++) {
        std::getline(file, line);
        std::stringstream linestream(line);
        for (int j = 0; j < 4; j++) {
            std::getline(linestream, value, ',');
        }
        std::getline(linestream, value, ',');
        grid[x][y][z] = std::stod(value);
        grid[x][y][z] = log10(grid[x][y][z]);
        x++;
        if (x >= 86) {
            y++;
            x = 0;
            if (y >= 71) {
                z++;
                y = 0;
            }
        }
    }

    file.close();
//    for (int i = 0; i < X; i++) {
//        for (int j = 0; j < Y; j++) {
//            for (int k = 0; k < Z; k++) {
//                std::cout << i << " " << j << " " << k << ": " << grid[i][j][k] << std::endl;
//            }
//        }
//    }

    MultilinearInterpolator<3> MLI(&grid[0][0][0], &minmax[0][0], &dims[0]);

    std::ofstream outfile;
    outfile.open("interp");

    double coord[3];
    coord[2] = -0.17;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (int i = 0; i < 100; i++) {
        for (int j = 0; j < 100; j++) {
//            std::cout << i << " " << j << std::endl;
            coord[0] = 2 + i * (8-2)/100.;
            coord[1] = -4 + j * 8 / 100.;
//            std::cout << coord[0] << " " << coord[1] << " " << coord[2] << std::endl;
            double result = MLI.interpolate(&coord[0]);

            outfile << coord[0] << " " << coord[1] << " " << result << std::endl;
            std::cout << std::endl;
        }
    }
    std::cout << "Done" << std::endl;
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time to complete = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

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
