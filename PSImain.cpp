//
// Created by Stefan Lüders on 26/07/2021.
//

#include "PSI.h"
#include "CoolConst.h"
#include "CoolSimplex.h"
#include <string>
#include <iostream>
#include <float.h>
#include <cstdlib>
#include <stdlib.h>
#include <cassert>
#include <math.h>
#include <time.h>


void test_knn() {
  
  srand (time(NULL));
  
  psi_init();
  psi_read_points("fulldata2.csv");
  psi_construct_btree();
  
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


void test_snn() {
  
  srand (time(NULL));
  
  psi_init();
  psi_read_points("fulldata2.csv");
  psi_construct_btree();
  
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


void test_projective_simplex(int adaptive) {
  
  srand (time(NULL));
  
  psi_init();
  psi_read_points("fulldata2.csv");
  psi_construct_btree();
  
  int count = 0;
  int iterations = 1000;
  int k = 30;
  double factor = 1.5;
  double max_steps = 5;
  
  
  // for testing projective simplex algorithm
  for (int i = 0; i < iterations; i++) {
    double target[3];
    target[0] = 2 + float(rand()) / RAND_MAX * 5;
    target[1] = -3 + float(rand()) / RAND_MAX * 5;
    target[2] = -2 + float(rand()) / RAND_MAX * 2;

//    target[0] = 5;
//    target[1] = 0;
//    target[2] = -1;
    
    int * simplex;
    if (adaptive) {
      simplex = psi_adaptive_projective_simplex_algorithm(target, k, factor, max_steps);
    } else {
      int * neighbours = psi_find_k_nearest_neighbor(target, k);
      simplex = psi_projective_simplex_algorithm(neighbours, target, k);
    }
//    std::cout << target[0] << " " << target[1] << " " << target[2] << ":\t";
//    for (int i = 0; i < k; i++) {
//      std::cout << neighbours[i] << " ";
//    }
//    std::cout << std::endl;
    
//    std::cout << i << ": ";
    if (simplex == nullptr) {
//      std::cout << "nullptr" << std::endl;
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
          
          sobj.points[j] = &points[j];
        }
        
        sobj.construct_T_inv();
        
        double * bary = sobj.convert_to_bary(target);
        int inside = sobj.check_bary(bary);
        
        if (!inside) {
          std::cout << i << "not inside" << std::endl;
        } else {
//          for (int j = 0; j < DIP_DIMS+1; j++ ) {
//            std::cout << bary[j] << " ";
//          }
//          std::cout << std::endl;
        }
      }
//      for (int j = 0; j < DIP_DIMS+1; j++) {
//        std::cout << simplex[j] << " ";
//      }
//      std::cout << std::endl;
    }
  }
  std::cout << "Loaded points: " << get_nlim() << std::endl;
  std::cout << "k = " << k << std::endl;
  if (adaptive) {
    std::cout << "factor: " << factor << std::endl;
    std::cout << "max_steps = " << max_steps << std::endl;
  }
  std::cout << iterations << " iterations" << std::endl;
  std::cout << "# of no solutions: " << count << std::endl;
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
  
  
  
  return 0;
}















