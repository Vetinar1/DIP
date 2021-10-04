//
// Created by Stefan Lüders on 26/07/2021.
//

#include "PSI.h"
#include "CoolConst.h"
#include "CoolCool.h"
#include "CoolSimplex.h"
#include <string>
#include <iostream>
#include <float.h>
#include <cstdlib>
#include <stdlib.h>
#include <cassert>
#include <math.h>
#include <time.h>
#include <chrono>
#include <fstream>
#include "CoolManager.h"

const int N_EVAL = 1000;


void read_eval(std::string cool_file, double ** coords, double ** vals, int apply_log_to_vals) {
  std::ifstream file;
  std::string line;
  std::string value;
  
  
  /* Read points */
  file.open(cool_file);
  if (!file.is_open()) {
    std::cerr << "Error reading " << cool_file << std::endl;
    abort();
  } else {
    std::cout << "Reading " << cool_file << " (max " << DIP_NMAX << " lines)" << std::endl;
  }

#ifdef DIP_POINTS_HEADER_SKIP
  std::getline(file, line);
#endif
  
  int n = 0;
  for (int i = 0; i < N_EVAL; i++) {
    std::getline(file, line);
    std::stringstream linestream(line);
    for (int j = 0; j < DIP_DIMS; j++) {     // DIP_DIMS coordinates
      std::getline(linestream, value, ',');
      coords[i][j] = std::stod(value);
    }
    for (int j = 0; j < DIP_VARNR; j++) {   // Htot, Ctot
      std::getline(linestream, value, ',');
      if (apply_log_to_vals) {
        vals[i][j] = log10(std::stod(value));
      } else {
        vals[i][j] = std::stod(value);}
    }
    
    n++;
    if (file.peek() == EOF) {
      break;
    }
  }
//  N_LIM = n;
  
  file.close();
}


//void test_psi() {
//  // for comparison
//
//  double ** evalcoords = new double * [N_EVAL];
//  double ** evalvals   = new double * [N_EVAL];
//
//  for (int i = 0; i < N_EVAL; i++) {
//    evalcoords[i] = new double[DIP_DIMS];
//    evalvals[i]   = new double[2]; // Htot, Ctot
//  }
//  std::string evalfiles[5] = {
//      "eval_01/tempdata.csv",
//      "eval_02/tempdata.csv",
//      "eval_03/tempdata.csv",
//      "eval_04/tempdata.csv",
//      "eval_05/tempdata.csv"
//  };
//
//  read_eval(evalfiles[DIP_DIMS-2], evalcoords, evalvals, 1);
//
//  srand (time(NULL));
//
////  double clamp_mins[DIP_DIMS];
////  double clamp_maxs[DIP_DIMS];
//  double clamp_mins[DIP_DIMS] = {4, -2, -1};//, -4, 7};//, 18.5};
//  double clamp_maxs[DIP_DIMS] = {6, 2, 0};//, 2, 9};//, 22.5};
//
//  std::string folder = "psi_13/";
//
//  std::cout << "Creating PSI object... " << std::endl;
//  PSI * psi = new PSI;
//
//  psi->set_clamp_values(clamp_mins, clamp_maxs);
//
//  std::cout << "Reading files... ";
//
//  std::chrono::steady_clock::time_point time1 = std::chrono::steady_clock::now();
//
//  psi->read_files(
//        folder + "z0.0.points"
////          folder + "tempdata.csv",
//  );
//
//  std::chrono::steady_clock::time_point time2 = std::chrono::steady_clock::now();
//  std::cout << "Time to read files = " << std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count() << "[ms]" << std::endl;
//
//  std::cout << "Done" << std::endl << "Beginning interpolation... " << std::endl;
//
//  const int iterations = 1000;
////  double qualities[iterations];
//
//  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
//
//  for (int i = 0; i < 1000; i++) {
//    double * target = evalcoords[i];
//    double * interps = psi->interpolate(target);
//  }
//  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
//  std::cout << "Done" << std::endl;
//  std::cout << "Time to complete = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
//
//}


void test_cm_psi() {
    
    double ** evalcoords = new double * [N_EVAL];
    double ** evalvals   = new double * [N_EVAL];
    
    for (int i = 0; i < N_EVAL; i++) {
        evalcoords[i] = new double[DIP_DIMS];
        evalvals[i]   = new double[DIP_VARNR];
    }
    
    read_eval("test_data/z0.42.points", evalcoords, evalvals, 0);
    
  std::cout << "Creating CoolManager object.." << std::endl;
  CoolManager * CM = new CoolManager(0.0, 0.5, "test_data/mapfile");

  std::cout << "Setting clamps..." << std::endl;
//  CM->set_clamp_values(clamp_mins, clamp_maxs);

  std::cout << "Interpolating 10000 points using CoolManager" << std::endl;
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  int count = 0;
  for (int i = 0; i < 922; i++) {
      std::cout << i << std::endl;
        double * coord = evalcoords[i];
        double * interp = CM->interpolate(coord, 0.42);
        
        if (interp != nullptr) {
            for (int j = 0; j < DIP_VARNR; j++) {
                std::cout << evalvals[i][j] << "\t";
            }
            std::cout << std::endl;
            for (int j = 0; j < DIP_VARNR; j++) {
                std::cout << interp[j] << "\t";
            }
            std::cout << std::endl;
            for (int j = 0; j < DIP_VARNR; j++) {
                std::cout << interp[j] / evalvals[i][j] << "\t";
            }
            std::cout << std::endl;
            std::cout << std::endl;
        }
        else {
            count++;
        }
  }
  
  std::cout << count << std::endl;
  std::cout << "Done" << std::endl;
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Time to complete = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    
}


int main() {
  // T: 1.3, 9.7
  // nH: -4.8, 4.8
  // Z: -2.3, 1.3
  
//  test_snn();
//  test_knn();
//  test_psi();
    test_cm_psi();
  
  return 0;
}















