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
#include "CoolMLI.h"

int main() {
    std::string mode = "block2d";
    std::cout << mode << std::endl;

//    if (mode == "slice3d") {
//        std::ifstream mapfile;
//        mapfile.open("../slice3d/mapfile");

//        if (!mapfile.is_open()) {
//            std::cerr << "Error reading mapfile" << std::endl;
//            return 1;
//        }

//        std::map<double, std::string> filenames;
//        double key;
//        std::string value;
//        while (mapfile >> key >> value) {
//            filenames[key] = "../slice3d/" + value;
//            std::cout << key << " " << value << " " << filenames[key] << std::endl;
//        }
//        mapfile.close();
//        CoolManager<1232, 3, 2444> CM(3.9, 4, filenames);

//        std::ofstream outfile;
//        outfile.open("interp");
//        double coord[2]; // T, nH
//        double z = 3.95;
//        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
//        // Parallel slice
// //        for (int i = 0; i < 100; i++) {
// //            for (int j = 0; j < 100; j++) {
// //                std::cout << i << " " << j << std::endl;
// //                coord[0] = 2 + i * (8-2)/100.;
// //                coord[1] = -4 + j * 8 / 100.;
// //                double interp = CM.interpolate(coord, z);
// //
// //                outfile << coord[0] << " " << coord[1] << " " << interp << std::endl;
// //                std::cout << std::endl;
// //            }
// //        }
//        // Orthogonal slice
//        coord[1] = 3;
//        z = 4;
//        for (int i = 0; i < 100; i++) {
//            for (int j = 0; j < 100; j++) {
//                coord[0] = 2 + j * (8-2)/100.;
//                z = 4 - i * 4 / 100.;
//                std::cout << coord[0] << " " << z << std::endl;
//                double interp = CM.interpolate(coord, z);

//                outfile << coord[0] << " " << z << " " << interp << std::endl;
//                std::cout << std::endl;
//            }
//        }
//        std::cout << "Done" << std::endl;
//        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
//        std::cout << "Time to complete = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
//        outfile.close();
//        return 0;
//    }

    if (mode == "block2d") {
        std::cout << "Initializing cool object... " << std::endl;
        //    Cool<35731, 4, 1002570> cool;
//        Cool<5994, 3, 38602> cool;
//        Cool<1050, 2, 2082> cool;
//        Cool<10000, 2, 10000> cool;
        Cool * cool = new Cool;
//        Cool<980, 2, 1940> cool;

        std::cout << "Reading files... ";
        cool->read_files("../gasoline_header2_tri/z0.0.points",
                    "../gasoline_header2_tri/z0.0.tris",
                    "../gasoline_header2_tri/z0.0.neighbors");
//        cool->read_files("../data2d/data.csv", "../data2d/dtri.csv", "../data2d/dneighbours.csv");
//        cool.read_files("../slice3d/z3.9.points", "../slice3d/z3.9.tris", "../slice3d/z3.9.neighbors");
//        cool.read_files("../slice3d/z3.9.points", "../slice3d/z3.9.tris", "../data3d/z3.9.neighbors");

        std::cout << "Done" << std::endl << "Constructing ball tree... " << std::flush;
        cool->construct_btree();
        // std::cout << "Done" << std::endl << "Saving ball tree... ";
        // cool->save_btree("tree");
        std::cout << "Done" << std::endl << "Beginning interpolation... " << std::endl;

        std::ofstream outfile;
        outfile.open("interp");
//        Parameter space:
//        T:	[2, 9]		Margins: 0.1
//        nH:	[-9, 4]		Margins: 0.1
//        Radiation background parameters:
//        SFR:	[-5, 3]		Margins: 0.1		Source file: spectra/SFR
//        old:	[6, 12]		Margins: 0.1		Source file: spectra/old
        double coord[D];
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                for (int k = 0; k < 10; k++) {
                    for (int l = 0; l < 10; l++) {
//                         std::cout << i << " " << j << std::endl;
                        coord[0] = 2 + i * (9-2)/10.;
//                        coord[1] = -9 + j * (4+9) / 10.;
                        coord[1] = -3 + j * (4+3) / 10.;
                        coord[2] = -5 + k * (8) / 10.;
                        coord[3] = 6 + (12-6) / 10.;

                        double interp = cool->interpolate(coord);

//                         outfile << coord[0] << " " << coord[1] << " " << interp << std::endl;
//                         std::cout << std::endl;

                    }
                }
            }
        }
        std::cout << "Done" << std::endl;
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Time to complete = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
        std::cout << "Avg = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / pow(100, D) << "[ms]" << std::endl;
        std::cout << "Avg flips: " << cool->avg_flips << std::endl;
        return 0;
    }

   if (mode == "multilinear") {
       int dims[4] = {71, 14, 7, 9};
       double * grid = new double[71*14*7*9];

       for (int i = 0; i < 71*14*7*9; i++) {
           grid[i] = 0;
       }

       // T, nH, old, SFR
       double minmax[4][2] = {
               {2, 9},
               {-9, 4},
               {6, 12},
               {-5, 3}
       };

       std::ifstream file;
       std::string line;
       std::string value;
       file.open("../gasoline_header2_grid/grid_gasoline_header2.csv");
       if (!file.is_open()) {
           std::cerr << "Error reading grid file" << std::endl;
           return 2;
       }
       for (int i = 0; i < 61760; i++) {
           std::getline(file, line);
           if (i == 0) {
               continue;
           }
           std::cout << i << std::endl;
           std::stringstream linestream(line);
           double coords[4];
           for (int j = 0; j < 4; j++) {
               std::getline(linestream, value, ',');
               std::cout << j << " " << value << " -> ";
               coords[j] = (std::stod(value) - minmax[j][0]) * (dims[j] - 1) / (minmax[j][1] - minmax[j][0]);
               std::cout << coords[j] << std::endl;
           }
           std::cout << 9 * (7 * ( 14 * (coords[0]) + coords[1]) + coords[2]) + coords[3] << std::endl;
           std::cout << std::endl;
           std::getline(linestream, value, ',');
           if (grid[(int) std::round(9 * (7 * ( 14 * (coords[0]) + coords[1]) + coords[2]) + coords[3])] != 0) {
               std::cerr << "reassignment" << std::endl;
               abort();
           }
           grid[
                   (int) std::round(9 * (7 * ( 14 * (coords[0]) + coords[1]) + coords[2]) + coords[3])
           ] = std::stod(value);
       }

       file.close();

       MultilinearInterpolator MLI(grid, &minmax[0][0], &dims[0]);

       std::ofstream outfile;
       outfile.open("interp");

       double coord[D];
//       coord[2] = -0.17;
       std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
//       for (int i = 0; i < 100; i++) {
//           for (int j = 0; j < 100; j++) {
//               std::cout << i << " " << j << std::endl;
//               coord[0] = 2 + i * (8-2)/100.;
//               coord[1] = -4 + j * 8 / 100.;
//   //            std::cout << coord[0] << " " << coord[1] << " " << coord[2] << std::endl;
//               double result = MLI.interpolate(&coord[0]);
//
//               outfile << coord[0] << " " << coord[1] << " " << result << std::endl;
//               std::cout << std::endl;
//           }
//       }
       for (int i = 0; i < 10; i++) {
           for (int j = 0; j < 10; j++) {
               for (int k = 0; k < 10; k++) {
                   for (int l = 0; l < 10; l++) {
//                         std::cout << i << " " << j << std::endl;
                       coord[0] = 2 + i * (9-2)/10.;
//                        coord[1] = -9 + j * (4+9) / 10.;
                       coord[1] = -3 + j * (4+3) / 10.;
                       coord[2] = -5 + k * (8) / 10.;
                       coord[3] = 6 + (12-6) / 10.;

                       double interp = MLI.interpolate(coord);

//                         outfile << coord[0] << " " << coord[1] << " " << interp << std::endl;
//                         std::cout << std::endl;

                   }
               }
           }
       }
       std::cout << "Done" << std::endl;
       std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
       std::cout << "Time to complete = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

       return 0;

   }
}
