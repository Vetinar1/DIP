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
void psi_set_clamp_values(double * cmins, double * cmaxs);
int psi_read_points(std::string cool_file, int apply_log_to_vals);
int psi_construct_btree(double ** points);
PSIBallTree * psi_construct_btree_recursive(double ** base, int * indices, int n);
int * psi_find_k_nearest_neighbor(double * target, int k);
int * psi_find_k_nearest_neighbor(double * target, double ** points, PSIBallTree * btree, int k);
int * psi_find_k_nearest_neighbour_bruteforce(double * target, int k);
double get_dist(double * target, int index);
double dot(double * v, double * w);
int * psi_projective_simplex_algorithm(int * neighbors, double * target, int k);
double get_coord(int i, int j);
double get_val(int i, int j);
int * psi_adaptive_projective_simplex_algorithm(double * target, int k, double factor, int max_steps);
int get_nlim();
double * psi_interpolate(double * target);
#endif //DIP_PSI_H
