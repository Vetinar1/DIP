//
// Created by Stefan Lueders on 15.09.20.
//

#include <iostream>
#include <chrono>
#include "cool.h"
#include "CoolManager.h"
#include <climits>
#include <bitset>
#include <math.h>
#include <string>
#include <map>
#include <iostream>

int main() {
    std::string mode = "block2d";

    if (mode == "slice3d") {
        std::ifstream mapfile;
        mapfile.open("../slice3d/mapfile");

        if (!mapfile.is_open()) {
            std::cerr << "Error reading mapfile" << std::endl;
            return 1;
        }

        std::map<double, std::string> filenames;
        double key;
        std::string value;
        while (mapfile >> key >> value) {
            filenames[key] = "../slice3d/" + value;
            std::cout << key << " " << value << " " << filenames[key] << std::endl;
        }
        mapfile.close();
        CoolManager<1293, 3, 2562> CM(3.9, 4, filenames);
    }

    if (mode == "block2d") {
        std::cout << "Initializing cool object... ";
        //    Cool<35731, 4, 1002570> cool;
//        Cool<5994, 3, 38602> cool;
        Cool<1050, 2, 2082> cool;
//        Cool<980, 2, 1940> cool;

        std::cout << "Done" << std::endl << "Reading files... ";
        cool.read_files("../data2d/data.csv", "../data2d/dtri.csv", "../data2d/dneighbours.csv");
//        cool.read_files("../slice3d/z3.9.points", "../slice3d/z3.9.tris", "../data3d/z3.9.neighbors");

        std::cout << "Done" << std::endl << "Constructing ball tree... " << std::flush;
        cool.construct_btree();
        std::cout << "Done" << std::endl << "Saving ball tree... ";
        cool.save_btree("tree");
        std::cout << "Done" << std::endl << "Beginning interpolation... " << std::endl;

        std::ofstream outfile;
        outfile.open("interp");

        double coord[2];
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        for (int i = 0; i < 100; i++) {
            for (int j = 0; j < 100; j++) {
                std::cout << i << " " << j << std::endl;
                coord[0] = 2 + i * (8-2)/100.;
                coord[1] = -4 + j * 8 / 100.;
                double interp = cool.interpolate(coord);

                outfile << coord[0] << " " << coord[1] << " " << interp << std::endl;
                std::cout << std::endl;
            }
        }
        std::cout << "Done" << std::endl;
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Time to complete = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
        std::cout << "Avg flips: " << cool.avg_flips << std::endl;
        return 0;
    }
//
//    if (mode == "multilinear") {
//        const int X = 86;
//        const int Y = 71;
//        const int Z = 10;
//        int dims[3] = {X, Y, Z};
//        double grid[X][Y][Z];
//
//        double minmax[3][2] = {
//                {1, 9.5},
//                {-8, 6},
//                {-4, 0.5}
//        };
//
//
//        std::ifstream file;
//        std::string line;
//        std::string value;
//        file.open("../repr1_schaye.csv");
//
//        int x = 0;
//        int y = 0;
//        int z = 0;
//        for (int i = 0; i < X*Y*Z; i++) {
//            std::getline(file, line);
//            std::stringstream linestream(line);
//            for (int j = 0; j < 4; j++) {
//                std::getline(linestream, value, ',');
//            }
//            std::getline(linestream, value, ',');
//            grid[x][y][z] = std::stod(value);
//            grid[x][y][z] = log10(grid[x][y][z]);
//            x++;
//            if (x >= 86) {
//                y++;
//                x = 0;
//                if (y >= 71) {
//                    z++;
//                    y = 0;
//                }
//            }
//        }
//
//        file.close();
//    //    for (int i = 0; i < X; i++) {
//    //        for (int j = 0; j < Y; j++) {
//    //            for (int k = 0; k < Z; k++) {
//    //                std::cout << i << " " << j << " " << k << ": " << grid[i][j][k] << std::endl;
//    //            }
//    //        }
//    //    }
//
//        MultilinearInterpolator<3> MLI(&grid[0][0][0], &minmax[0][0], &dims[0]);
//
//        std::ofstream outfile;
//        outfile.open("interp");
//
//        double coord[3];
//        coord[2] = -0.17;
//        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
//        for (int i = 0; i < 100; i++) {
//            for (int j = 0; j < 100; j++) {
//                std::cout << i << " " << j << std::endl;
//                coord[0] = 2 + i * (8-2)/100.;
//                coord[1] = -4 + j * 8 / 100.;
//    //            std::cout << coord[0] << " " << coord[1] << " " << coord[2] << std::endl;
//                double result = MLI.interpolate(&coord[0]);
//
//                outfile << coord[0] << " " << coord[1] << " " << result << std::endl;
//                std::cout << std::endl;
//            }
//        }
//        std::cout << "Done" << std::endl;
//        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
//        std::cout << "Time to complete = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
//
//        return 0;
//
//    }
}
