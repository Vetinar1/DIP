//
// Created by Stefan Lüders on 26/07/2021.
//

#ifndef DIP_PSI_H
#define DIP_PSI_H

#include <string>
#include "PSIBallTree.h"
#include <queue>
#include <iostream>
#include "Const.h"

struct distpoint {
    int index;
    double dist;
    bool operator<(const distpoint& rhs) const {
        return dist < rhs.dist;
    }
};


class PSI {
private:
    double ** coords;
    double ** vals;
    double CLAMP_MIN[DIP_DIMS];
    double CLAMP_MAX[DIP_DIMS];
    int N_LIM;
    PSIBallTree * btree;
    
    double dot(double * v, double * w);
    PSIBallTree * construct_btree_recursive(double ** base, int * indices, int n);
    void find_k_nearest_neighbor_recursive(PSIBallTree * root, const double * target, double ** base,
                                            std::priority_queue <distpoint> * Q, int k);
    
public:
    PSI() {
      coords = new double * [DIP_NMAX];
      vals = new double * [DIP_NMAX];
  
      for (int i = 0; i < DIP_NMAX; i++) {
        coords[i] = new double[DIP_DIMS];
        vals[i]   = new double[DIP_VARNR];
      }
    };
    
    int read_files(std::string cool_file);
    int get_nlim();
    void reset();
    double get_coord(int i, int j);
    double get_val(int i, int j);
    double get_dist(double * target, int index);
    void set_clamp_values(double * cmins, double * cmaxs);
    int construct_btree(double ** points);
    int * find_k_nearest_neighbor(double * target, double ** points, PSIBallTree * btree, int k);
    int * find_k_nearest_neighbor(double * target, int k);
    int * projective_simplex_algorithm(int * neighbors, double * target, int k, int smart_nn);
    int * adaptive_projective_simplex_algorithm(double * target, int k, double factor, int max_steps, int smart_nn);
    double * interpolate(double * target, int smart_fallback);
};

#endif //DIP_PSI_H
