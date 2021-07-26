//
// Created by Stefan Lüders on 26/07/2021.
//

#ifndef DIP_PSI_H
#define DIP_PSI_H

#include <string>
#include "PSIBallTree.h"

void psi_init();
int psi_read_points(std::string cool_file);
int psi_construct_btree();
PSIBallTree * psi_construct_btree_recursive(int * indices, int n);
PSIBallTree * psi_find_nearest_neighbour_recursive(PSIBallTree * root, const double * target,
                                                   PSIBallTree * best, double * min_dist2);
void psi_find_nearest_neighbour_recursive(PSIBallTree * root, const double * target, double * min_dist2);
int psi_find_nearest_neighbour(double * target);
int * psi_find_nearest_neighbour(double * target, int k);
int psi_find_nearest_neighbour_bruteforce(double * target);
double get_dist(double * target, int index);
#endif //DIP_PSI_H
