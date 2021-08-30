//
// Created by Stefan Lüders on 26/07/2021.
//

#ifndef DIP_PSI_H
#define DIP_PSI_H

#include <string>
#include "PSIBallTree.h"
#include <queue>
#include <iostream>

struct distpoint {
    int index;
    double dist;
    bool operator<(const distpoint& rhs) const {
      return dist < rhs.dist;
    }
};

void psi_init();
int psi_read_points(std::string cool_file);
int psi_construct_btree();
PSIBallTree * psi_construct_btree_recursive(int * indices, int n);
PSIBallTree * psi_find_nearest_neighbour_recursive(PSIBallTree * root, const double * target,
                                                   PSIBallTree * best, double * min_dist2);
void psi_find_nearest_neighbour_recursive(PSIBallTree * root, const double * target,
                                          std::priority_queue<distpoint> * Q, int k);
int psi_find_nearest_neighbour(double * target);
int psi_find_nearest_neighbour_bruteforce(double * target);
int * psi_find_k_nearest_neighbor(double * target, int k);
int * psi_find_k_nearest_neighbour_bruteforce(double * target, int k);
double get_dist(double * target, int index);
double dot(double * v, double * w);
int * psi_projective_simplex_algorithm(int * neighbors, double * target, int k);
double get_coord(int i, int j);
int * psi_adaptive_projective_simplex_algorithm(double * target, int k, double factor, int max_steps);
int get_nlim();
double get_average_executions();
#endif //DIP_PSI_H
