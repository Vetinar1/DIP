//
// Created by Stefan Lüders on 26/07/2021.
//

#include "PSI.h"
#include "CoolConst.h"
#include <string>
#include <iostream>
#include <float.h>
#include <cstdlib>
#include <stdlib.h>
#include <cassert>
#include <math.h>
#include <time.h>

int main() {
  // T: 1.3, 9.7
  // nH: -4.8, 4.8
  // Z: -2.3, 1.3
  
  srand (time(NULL));
  
  psi_init();
  psi_read_points("fulldata2.csv");
  psi_construct_btree();
  
  int count = 0;
  
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
  return 0;
}