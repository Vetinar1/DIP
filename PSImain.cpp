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

const int N_EVAL = 1000;


void test_knn() {
  
  srand (time(NULL));
  
  psi_init();
  psi_read_points("fulldata2.csv", 0);
  
  int count = 0;
  
  
  //  // for testing k nearest neighbor code
  for (int i = 0; i < 100; i++) {
    double target[3];
    target[0] = 2 + float(rand()) / RAND_MAX * 5;
    target[1] = -3 + float(rand()) / RAND_MAX * 5;
    target[2] = -2 + float(rand()) / RAND_MAX * 2;
    
    std::cout << std::endl;
    std::cout << i << " " << target[0] << " " << target[1] << " " << target[2] << std::endl;
    
    int k = 50;
    int * knntree = psi_find_k_nearest_neighbor(target, k);
    int * knnbf = psi_find_k_nearest_neighbour_bruteforce(target, k);
    
    int different = 0;
    for (int j = 0; j < k; j++) {
      if (knntree[j] != knnbf[j]) {
        different = 1;
        break;
      }
    }
    if (different) {
      count += 1;
      for (int j = 0; j < k; j++) {
        std::cout << i << " neighbor " << j << ": " << knntree[j] << " / " << knnbf[j] << std::endl;
      }
    }
  }
  std::cout << count << std::endl;
}


#if 0
void test_snn() {
  
  srand (time(NULL));
  
  psi_init();
  psi_read_points("fulldata2.csv");
  
  int count = 0;
  
  // for testing single nearest neighbor code
  for (int i = 0; i < 1000; i++) {
    double target[3];
    target[0] = 2 + float(rand()) / RAND_MAX * 5;
    target[1] = -3 + float(rand()) / RAND_MAX * 5;
    target[2] = -2 + float(rand()) / RAND_MAX * 2;
    
    int indextree = psi_find_nearest_neighbour(target);
    int indexbf = psi_find_nearest_neighbour_bruteforce(target);
    
    if (indextree != indexbf) {
      count += 1;
      double dist1 = get_dist(target, indextree);
      double dist2 = get_dist(target, indexbf);
      std::cout << i << ": " << indextree << " " << dist1 << " / " << indexbf << " " << dist2 << std::endl;
    }
  }
  std::cout << count << std::endl;
}
#endif


void test_quality() {
  // Verify that the regular n-simplex has quality 1
  Simplex simplex;
  
  // Construct regular simplex out of basis vector + 1 additional vector, see wikipedia
  for (int i = 0; i < DIP_DIMS; i++) {
    simplex.points[i] = new Point;
    
    for (int j = 0; j < DIP_DIMS; j++) {
      simplex.points[i]->coords[j] = 0;
      
      if (i == j) {
        simplex.points[i]->coords[j] = 1 / sqrt(2);
      }
    }
  }
  
  // Remaining point
  simplex.points[DIP_DIMS] = new Point;
  for (int i = 0; i < DIP_DIMS; i++) {
    simplex.points[DIP_DIMS]->coords[i] = (1 + sqrt(DIP_DIMS+1)) / (DIP_DIMS * sqrt(2));
  }
  
  simplex.construct_T_inv();
  double quality = simplex.get_quality();
  std::cout << "Simplex has quality: " << quality << std::endl;
}


void test_dip() {
  // for comparison
  
  srand (time(NULL));
  
//  double clamp_mins[DIP_DIMS];
//  double clamp_maxs[DIP_DIMS];
  double clamp_mins[DIP_DIMS] = {4, -2};//, -1, -4, 7};//, 18.5};
  double clamp_maxs[DIP_DIMS] = {6, 2};//, 0, 2, 9};//, 22.5};
  
  std::cout << "Creating Cool object... " << std::endl;
  Cool * cool = new Cool;
  
  cool->set_clamp_values(clamp_mins, clamp_maxs);
  
  std::cout << "Reading files... ";
  
  std::chrono::steady_clock::time_point time1 = std::chrono::steady_clock::now();
  
  cool->read_files(
//      "synthetic/synthetic_5d.csv",
//      "synthetic/synthetic_5d.tris",
//      "synthetic/synthetic_5d.neighbors"
        "psi_02/z0.0.points",
        "psi_02/z0.0.tris",
        "psi_02/z0.0.neighbors"
//          "psi_11/tempdata.csv",
//          "psi_11/tempdata.tris",
//          "psi_11/tempdata.neighbors"
  );
  
  std::chrono::steady_clock::time_point time2 = std::chrono::steady_clock::now();
  std::cout << "Time to read files = " << std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count() << "[ms]" << std::endl;
  
  std::cout << "Done" << std::endl << "Constructing ball tree... " << std::flush;
  cool->construct_btree();
  
  std::cout << "Done" << std::endl << "Beginning interpolation... " << std::endl;
  
  const int iterations = 1000;
  double interps[iterations];
  
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  std::cout << "Time to build tree = " << std::chrono::duration_cast<std::chrono::milliseconds>(begin - time2).count() << "[ms]" << std::endl;
  
  for (int i = 0; i < 1000; i++) {
    double target[DIP_DIMS];
    target[0] = 4 + float(rand()) / RAND_MAX * 2;
    target[1] = -2 + float(rand()) / RAND_MAX * 4;
    target[2] = -1 + float(rand()) / RAND_MAX * 1;
//    target[3] = -4 + float(rand()) / RAND_MAX * 6;
//    target[4] = 7 + float(rand()) / RAND_MAX * 2;
//    target[5] = 18.5 + float(rand()) / RAND_MAX * 4;
    interps[i] = cool->interpolate(target);
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Done" << std::endl;
  std::cout << "Time to complete = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
  std::cout << "Avg = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / pow(100, D) << "[ms]" << std::endl;
  std::cout << "Avg flips: " << cool->avg_flips << std::endl;
  std::cout << "Avg quality (+- stdev): " << cool->quality_avg << "+-" << cool->quality_stdev << std::endl;
  std::cout << interps[0] << std::endl;
}


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
    for (int j = 0; j < 2; j++) {   // Htot, Ctot
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


void test_projective_simplex(int adaptive) {
  
  double ** evalcoords = new double * [N_EVAL];
  double ** evalvals   = new double * [N_EVAL];
  
  for (int i = 0; i < N_EVAL; i++) {
    evalcoords[i] = new double[DIP_DIMS];
    evalvals[i]   = new double[2]; // Htot, Ctot
  }
  std::string evalfiles[5] = {
      "eval_01/tempdata.csv",
      "eval_02/tempdata.csv",
      "eval_03/tempdata.csv",
      "eval_04/tempdata.csv",
      "eval_05/tempdata.csv"
  };
  
  read_eval(evalfiles[DIP_DIMS-2], evalcoords, evalvals, 1);
  
  srand (time(NULL));
  
  psi_init();
  
  std::string folder = "psi_01/";
  std::string fname = "z0.0.points";
//  psi_read_points("psi_13/tempdata.csv");
  psi_read_points(folder.append(fname), 0);
//  psi_read_points("synthetic/synthetic_6d.csv");
  
  int count = 0;
  int count2 = 0;
  const int iterations = 1000;
  int k = 20;
  double factor = 2;
  double max_repetitions = 4;
  
  double qualities[iterations];
  double interpolations[iterations];
  
  
  // for testing projective simplex algorithm
  
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  for (int i = 0; i < iterations; i++) {
//    double target[DIP_DIMS];
    
    double * target = evalcoords[i];
//    target[0] = 4 + float(rand()) / RAND_MAX * 2;
//    target[1] = -2 + float(rand()) / RAND_MAX * 4;
//    target[2] = -1 + float(rand()) / RAND_MAX * 1;
//    target[3] = -4 + float(rand()) / RAND_MAX * 6;
//    target[4] = 7 + float(rand()) / RAND_MAX * 2;
//    target[5] = 18.5 + float(rand()) / RAND_MAX * 4;
    
    int * simplex;
    if (adaptive) {
      simplex = psi_adaptive_projective_simplex_algorithm(target, k, factor, max_repetitions);
    } else {
      int * neighbours = psi_find_k_nearest_neighbor(target, k);
      simplex = psi_projective_simplex_algorithm(neighbours, target, k);
    }
    
    if (simplex == nullptr) {
      count++;
    } else {
      
      // Verify
      if (simplex != nullptr) {
        Point points[DIP_DIMS+1];
        Simplex sobj;
        
        for (int j = 0; j < DIP_DIMS+1; j++) {
          for (int k = 0; k < DIP_DIMS; k++) {
            points[j].coords[k] = get_coord(simplex[j], k);
          }
          
          points[j].value = get_val(simplex[j], 0);
          sobj.points[j] = &points[j];
        }
        
        sobj.construct_T_inv();
        
        double * bary = sobj.convert_to_bary(target);
        int inside = sobj.check_bary(bary);
        
        if (!inside) {
          std::cout << i << "not inside" << std::endl;
          count2++;
        }
        
        qualities[i] = sobj.get_quality();
        interpolations[i] = 0;
        for (int j = 0; j < DIP_DIMS+1; j++) {
//          std::cout << interpolations[i] << "+=" << bary[j] << " * " << sobj.points[j]->value << std::endl;
          interpolations[i] += bary[j] * sobj.points[j]->value;
        }
      }
    }
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Loaded " << get_nlim() << " points in " << DIP_DIMS << " dimensions" << std::endl;;
  std::cout << "Built " << iterations << " simplices in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
  std::cout << "k = " << k << std::endl;
  if (adaptive) {
    std::cout << "adaptive factor: " << factor << std::endl;
    std::cout << "max_repetitions = " << max_repetitions << std::endl;
    std::cout << "average executions: " << get_average_executions() << std::endl;
  }
  std::cout << "# of no solutions: " << count << std::endl;
  std::cout << "# of wrong solutions: " << count2 << std::endl;
  
  double q_mean = 0;
  double i_mean = 0;
  double q_M2 = 0;
  double i_M2 = 0;
  double q_stdev = 0;
  double i_stdev = 0;
  double q_min = DBL_MAX;
  double i_min = DBL_MAX;
  double q_max = 0;
  double i_max = 0;
  std::ofstream q_file;
  q_file.open(folder.append("psi_qualities.csv"));
  q_file << "simplex quality,interpolation,actual value\n";
  
  for (int i = 0; i < iterations; i++) {
//    std::cout << evalvals[i][1] << " " << interpolations[i] << std::endl;
    double i_diff_abs = fabs(evalvals[i][1] - interpolations[i]);
    q_file << qualities[i] << "," << interpolations[i] << "," << evalvals[i][1] << "\n";
    
    if (qualities[i] < q_min) {
      q_min = qualities[i];
    } else if (qualities[i] > q_max) {
      q_max = qualities[i];
    }
    double q_delta = qualities[i] - q_mean;
    double i_delta = i_diff_abs - i_mean;
    q_mean += q_delta/(i+1);
    i_mean += i_delta/(i+1);
    q_M2 += q_delta * (qualities[i] - q_mean);
    i_M2 += i_delta * (i_diff_abs - i_mean);
  }
  q_stdev = q_M2 / (iterations - 1);
  i_stdev = i_M2 / (iterations - 1);
  std::cout << "Simplex quality: " << q_mean << "+-" << q_stdev << " (in [" << q_min << ", " << q_max << "])" << std::endl;
  std::cout << "Avg. abs. error: " << i_mean << "+-" << i_stdev << std::endl;
  q_file.close();
}


void test_projective_simplex() {
  test_projective_simplex(0);
}


int main() {
  // T: 1.3, 9.7
  // nH: -4.8, 4.8
  // Z: -2.3, 1.3
  
//  test_snn();
//  test_knn();
  test_projective_simplex(1);
  
//  test_dip();
  
  
  return 0;
}















