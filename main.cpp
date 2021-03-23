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

#include <iomanip>
#include <sstream>

// TODO Important: Dont go all the way to core edges. -> add to coolmanager clamps

int main() {
    std::string mode = "block";
//    std::cout << mode << std::endl;
//    std::cout << std::setprecision( std::numeric_limits<double>::digits10+2);

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

    if (mode == "block") {
        std::cout << "Initializing cool object... " << std::endl;
        //    Cool<35731, 4, 1002570> cool;
//        Cool<5994, 3, 38602> cool;
//        Cool<1050, 2, 2082> cool;
//        Cool<10000, 2, 10000> cool;
        Cool * cool = new Cool;
//        Cool<980, 2, 1940> cool;

        double clamp_mins[D] = {2, -9, -3};
        double clamp_maxs[D] = {9, 4, 1};
//        double clamp_mins[D] = {2.1, -8.9, -1.9};
//        double clamp_maxs[D] = {8.9, 3.9, -0.1};
        cool->set_clamp_values(clamp_mins, clamp_maxs);

        std::cout << "Reading files... ";


        cool->read_files(
                "../run45_gadget/z0.0.points",
                "../run45_gadget/z0.0.tris",
                "../run45_gadget/z0.0.neighbors"
//                "../complexity/mesh4/experiment.points",
//                "../complexity/mesh4/experiment.tris",
//                "../complexity/mesh4/experiment.neighbors"
        );
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
        for (int i = 0; i < 100; i++) {
            coord[0] = 2 + i * 7. / 100;
            coord[1] = 0;
            coord[2] = -2;

            std::cout << i << " " << coord[0] << " ";
            std::cout << cool->interpolate(coord);
            std::cout << std::endl;
//            std::cout << std::endl;


//            for (int j = 0; j < 10; j++) {
//                for (int k = 0; k < 10; k++) {
////                    for (int l = 0; l < 10; l++) {
////                        std::cout << i << " " << j << " " << k << " " << l << std::endl;
//                        coord[0] = 2 + i * (9-2)/100.;
////                        coord[1] = -9 + j * (4+9) / 10.;
//                        coord[1] = -3 + j * (4+3) / 100.;
////                        coord[2] = -5 + k * (8) / 10.;
//                        coord[2] = -1.9 + k * (1.8) / 10.;
////                        coord[3] = 6.1 + l * (11.9-6.1) / 10.;
////                        std::cout << coord[0] << " " << coord[1] << " " << coord[2] << std::endl;
//                        double interp = cool->interpolate(coord);
//
////                         outfile << coord[0] << " " << coord[1] << " " << interp << std::endl;
////                         std::cout << std::endl;
//
////                    }
//                }
//            }
        }
        std::cout << "Done" << std::endl;
        std::cout << "Encountered " << cool->nullpointers_encountered << " Nullpointers in Simplex traversal while " <<
                    " interpolating " << cool->interpolate_calls << " points" << std::endl;
                  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
//        std::cout << "Time to complete = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[us]" << std::endl;
        std::cout << "Time to complete = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
        std::cout << "Avg = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / pow(100, D) << "[ms]" << std::endl;
        std::cout << "Avg flips: " << cool->avg_flips << std::endl;
        return 0;
    }

   if (mode == "multilinear") {
       int dims[D] = {71, 14};
       int n_points = 1;
       for (int i = 0; i < D; i++) {
           n_points *= dims[i];
       }

       double * grid = new double[n_points];

       for (int i = 0; i < n_points; i++) {
           grid[i] = 0;
       }

       // T, nH, old, SFR
       double minmax[D][2] = {
               {2, 9},
               {-9, 4}
       };

       std::ifstream file;
       std::string line;
       std::string value;
       file.open("../complexity/grids/grid2.csv");
       if (!file.is_open()) {
           std::cerr << "Error reading grid file" << std::endl;
           return 2;
       }
       for (int i = 0; i < 995; i++) {
           std::getline(file, line);
           if (i == 0) {
               continue;
           }
           std::stringstream linestream(line);

           // get point coords (first D columns)
           double coords[D];
           for (int j = 0; j < D; j++) {
               std::getline(linestream, value, ',');
               coords[j] = (std::stod(value) - minmax[j][0]) * (dims[j] - 1) / (minmax[j][1] - minmax[j][0]);
           }

           // get point value (last column)
           std::getline(linestream, value, ',');

           // determine grid index
           double grid_index = coords[0];
           for (int j = 1; j < D; j++) {
               grid_index = dims[j] * grid_index + coords[j];
           }
//           if (grid[(int) std::round(9 * (7 * ( 14 * (coords[0]) + coords[1]) + coords[2]) + coords[3])] != 0) {
           if (grid[(int) std::round(grid_index)] != 0) {
               std::cerr << "reassignment" << std::endl;
               abort();
           }
           grid[
                   (int) std::round(grid_index)
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
       for (int i = 0; i < 100; i++) {
           for (int j = 0; j < 100; j++) {
//               for (int k = 0; k < 10; k++) {
//                   for (int l = 0; l < 10; l++) {
//                         std::cout << i << " " << j << std::endl;
                       coord[0] = 2 + i * (9-2)/100.;
//                        coord[1] = -9 + j * (4+9) / 10.;
                       coord[1] = -3 + j * (4+3) / 100.;
//                       coord[2] = -5 + k * (8) / 10.;
//                       coord[3] = 6 + (12-6) / 10.;

                       double interp = MLI.interpolate(coord);

//                         outfile << coord[0] << " " << coord[1] << " " << interp << std::endl;
//                         std::cout << std::endl;

//                   }
//               }
           }
       }
       std::cout << "Done" << std::endl;
       std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
       std::cout << "Time to complete = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[us]" << std::endl;

       return 0;

   }
}
