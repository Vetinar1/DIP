//
// Created by Stefan Lüders on 26/07/2021.
//

#include "CoolConst.h"
#include "PSIBallTree.h"
#include <algorithm>
#include <math.h>
#include <cassert>
#include <string>
#include <float.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <queue>
#include <vector>

static double ** coords;
static double ** vals;

static int N_LIM;

static PSIBallTree * btree; // root node

void psi_init() {
  coords = new double * [DIP_NMAX];
  vals   = new double * [DIP_NMAX];
  
  for (int i = 0; i < DIP_NMAX; i++) {
    coords[i] = new double[DIP_DIMS];
    vals[i]   = new double[DIP_DIMS];
  }
}


int psi_read_points(std::string cool_file) {
  std::ifstream file;
  std::string line;
  std::string value;
  
  
  /* Read points */
  file.open(cool_file);
  if (!file.is_open()) {
    std::cerr << "Error reading " << cool_file << std::endl;
    return 1;
  } else {
    std::cout << "Reading " << cool_file << std::endl;
  }

#ifdef DIP_POINTS_HEADER_SKIP
  std::getline(file, line);
#endif
  
  int n = 0;
  for (int i = 0; i < DIP_NMAX; i++) {
    std::getline(file, line);
    std::stringstream linestream(line);
    for (int j = 0; j < DIP_DIMS; j++) {     // DIP_DIMS coordinates
      std::getline(linestream, value, ',');
      coords[i][j] = std::stod(value);
    }
    for (int j = 0; j < DIP_VARNR; j++) {   // DIP_VARNR values
      std::getline(linestream, value, ',');
      vals[i][j] = std::stod(value);
    }
    
    n++;
    if (file.peek() == EOF) {
      break;
    }
  }
  N_LIM = n;
  
  file.close();
  
  return 0;
}


PSIBallTree * psi_construct_btree_recursive(int * indices, int n) {
  // Input validation + what if theres only one element left?
  if (n < 0) {
    std::cerr << "Illegal argument in construct_btree: n = " << n << std::endl;
  } else if (n == 1) {
    PSIBallTree * leaf = new PSIBallTree;
    leaf->pivot  = indices[0];
    leaf->lchild = nullptr;
    leaf->rchild = nullptr;
    leaf->radius = 0;
    
    return leaf;
  }
  
  // More than one element
  // 1. Find dimension of greatest spread and its median
  double largest_spread = -1;
  int lspread_dim = -1;
  double median = 0;
  
  for (int i = 0; i < DIP_DIMS; i++) {   // every dimension
    double * coord_i = new double[n];
    
    for (int j = 0; j < n; j++) {
      coord_i[j] = coords[indices[j]][i];
    }
    
    // sort values to determine spread and median
    std::sort(coord_i, coord_i + n);
    double dim_spread = coord_i[n-1] - coord_i[0];

#ifdef DIP_BALLTREE_CHECKS
    assert(dim_spread >= 0);
#endif
    
    if (dim_spread < largest_spread) {
      delete[] coord_i;
      continue;
    } else {
      largest_spread = dim_spread;
      lspread_dim = i;
    }
    
    if (n % 2 == 0) {
      median = coord_i[n/2];
    } else {
      median = (coord_i[(int) floor(n/2)] + coord_i[(int) ceil(n/2)]) / 2;
    }
    delete[] coord_i;
  }

#ifdef DIP_BALLTREE_CHECKS
  assert(largest_spread != -1 && lspread_dim != -1);
#endif
  
  // 2. Find pivot - point closest to median                       and
  // 3. Group points into sets to the left and right of average (L, R)
  
  PSIBallTree * node = new PSIBallTree;
  int * L = new int[n];
  int * R = new int[n];
  int lcount = 0;
  int rcount = 0;
  
  double min_dist = DBL_MAX;
  node->pivot = -1;
  
  for (int i = 0; i < n; i++) {
    double val = coords[indices[i]][lspread_dim];
    double dist = fabs(median - val);
    
    if (dist < min_dist) {
      // Add old pivot to L or R
      if (node->pivot != -1) {
        if (coords[node->pivot][lspread_dim] <= median) {
          L[lcount] = node->pivot;
          lcount++;
        } else {
          R[rcount] = node->pivot;
          rcount++;
        }
      }
      
      // Save new pivot
      min_dist = dist;
      node->pivot = indices[i];
      
    } else {
      if (val <= median) {
        L[lcount] = indices[i];
        lcount++;
      } else {
        R[rcount] = indices[i];
        rcount++;
      }
    }
  }

#ifdef DIP_BALLTREE_CHECKS
  assert(node->pivot != -1);
  assert(rcount != 0 || lcount != 0);
  assert(min_dist < DBL_MAX);
  for (int i = 0; i < lcount; i++) {
      assert(L[i] != node->pivot);
  }
  for (int i = 0; i < rcount; i++) {
      assert(R[i] != node->pivot);
  }
#endif
  
  // 4. Recurse on L and R
  if (lcount > 0) {
    node->lchild = psi_construct_btree_recursive(L, lcount);
  } else {
    node->lchild = nullptr;
  }
  if (rcount > 0) {
    node->rchild = psi_construct_btree_recursive(R, rcount);
  } else {
    node->rchild = nullptr;
  }

#ifdef DIP_BALLTREE_CHECKS
  // This should only happen if there is only one element left, and in that case this line shouldnt be reached
    assert(node->lchild != nullptr || node->rchild != node->rchild);
#endif
  
  // 5. Determine ball radius
  // squared distance to farthest element in subtree
  node->radius = 0;
  for (int i = 0; i < n; i++) {   // n elements of input array
    double dist = 0;
    // distance calculation
    for (int j = 0; j < DIP_DIMS; j++) {   // D coordinates
      dist += pow(coords[indices[i]][j] - coords[node->pivot][j], 2);
    }
    dist = sqrt(dist);
    
    if (dist > node->radius) {
      node->radius = dist;
    }
  }
  
  // 6. Bonus: Determine which child is closer
  if (node->lchild != nullptr && node->rchild != nullptr) {
    
    double lchild_dist2 = 0;
    double rchild_dist2 = 0;
    for (int i = 0; i < DIP_DIMS; i++) {
      lchild_dist2 += pow(coords[node->pivot][i] - coords[node->lchild->pivot][i], 2);
      rchild_dist2 += pow(coords[node->pivot][i] - coords[node->rchild->pivot][i], 2);
    }
  
    if (lchild_dist2 < rchild_dist2) {
      node->closechild = node->lchild;
      node->farchild = node->rchild;
    } else {
      node->closechild = node->rchild;
      node->farchild = node->lchild;
    }
    
  } else {
    // If one or both are nullpointers order doesnt matter
    node->closechild = node->lchild;
    node->farchild = node->rchild;
  }
  
  delete[] L;
  delete[] R;
  
  return node;
}


int psi_construct_btree() {
  // Construct array of indices
  int * indices = new int[N_LIM];
  for (int i = 0; i < N_LIM; i++) {
    indices[i] = i;
  }
  
  btree = psi_construct_btree_recursive(indices, N_LIM);
  delete[] indices;
  
  return 0;
}


PSIBallTree * psi_find_nearest_neighbour_recursive(PSIBallTree * root, const double * target,
                                                   PSIBallTree * best, double * min_dist2) {
  // Distance to current node
  double dist2 = 0;
  for (int i = 0; i < DIP_DIMS; i++) {
    dist2 += pow(target[i] - coords[root->pivot][i], 2);
  }
  dist2 = sqrt(dist2);
  
  // Recursion exit condition - target point is further outside the current node's ball than the distance
  // to the current closest neighbor -> there can't be a closer neighbor in this ball/subtree
  if (dist2 - root->radius >= sqrt(*min_dist2)) {
    return best;
  }
  
  // Perhaps the current point is the new closest?
  if (dist2 < *min_dist2) {
    *min_dist2 = dist2;
    best = root;
  }
  
  // Find closest child
  double ldist2 = 0;
  double rdist2 = 0;
  
  if (root->lchild != nullptr) {
    for (int i = 0; i < DIP_DIMS; i++) {
      ldist2 += pow(target[i] - coords[root->lchild->pivot][i], 2);
    }
  }
  if (root->rchild != nullptr) {
    for (int i = 0; i < DIP_DIMS; i++) {
      rdist2 += pow(target[i] - coords[root->rchild->pivot][i], 2);
    }
  }
  
  // Recurse into closest child first
  if (ldist2 <= rdist2) {
    if (root->lchild != nullptr) {
      best = psi_find_nearest_neighbour_recursive(root->lchild, target, best, min_dist2);
    }
    if (root->rchild != nullptr) {
      best = psi_find_nearest_neighbour_recursive(root->rchild, target, best, min_dist2);
    }
  } else {
    if (root->rchild != nullptr) {
      best = psi_find_nearest_neighbour_recursive(root->rchild, target, best, min_dist2);
    }
    if (root->lchild != nullptr) {
      best = psi_find_nearest_neighbour_recursive(root->lchild, target, best, min_dist2);
    }
  }
  
  return best;
}


void psi_find_nearest_neighbour_recursive(PSIBallTree * root, const double * target,
                                          std::priority_queue<int> * Q, int k) {
  /**
   * This function may look a bit odd because it was originally only supposed to find *one* nearest neighbour.
   */
  // Distance to current node and distance to farthest node in Q
  double root_dist = 0;
  double Q_dist = 0;
  for (int i = 0; i < DIP_DIMS; i++) {
    root_dist += pow(target[i] - coords[root->pivot][i], 2);
    Q_dist += pow(target[i] - coords[Q->top()][i], 2);
  }
  root_dist = sqrt(root_dist);
  Q_dist = sqrt(Q_dist);

  // Recursion exit condition - target point is further outside the current node's ball than the distance
  // to the current farthest nearest neighbor -> there can't be a closer neighbor in this ball/subtree
  if (root_dist - root->radius >= Q_dist) {
    return;
  }

  // The current node is closer than at least one node from the queue. Add it, keep queue length <= k
  Q->push(root->pivot);
  if (Q->size() > k) {
    Q->pop();
  }
  
  if (root->closechild != nullptr) {
    psi_find_nearest_neighbour_recursive(root->closechild, target, Q, k);
  }
  if (root->farchild != nullptr) {
    psi_find_nearest_neighbour_recursive(root->farchild, target, Q, k);
  }
}


int psi_find_nearest_neighbour(double * target) {
  double min_dist2 = DBL_MAX;
  PSIBallTree * best = nullptr;
  PSIBallTree * node = psi_find_nearest_neighbour_recursive(btree, target, best, &min_dist2);
  
  return node->pivot;
}


int * psi_find_nearest_neighbour(double * target, int k) {
  std::priority_queue<int> * Q;
  double min_dist2 = DBL_MAX;
  PSIBallTree * best = nullptr;
  psi_find_nearest_neighbour_recursive(btree, target, Q, k);
  
  assert(Q->size() == k);
  
  int * out = new int[k];
  for (int i = k-1; i > 0; i--) {
    out[i] = Q->top();
    Q->pop();
  }
  
  assert(Q->size() == 0);

  return out;
}


int psi_find_nearest_neighbour_bruteforce(double * target) {
  double min_dist2 = DBL_MAX;
  int min_index;
  
  for (int i = 0; i < N_LIM; i++) {
    double dist2 = 0;
    for (int j = 0; j < DIP_DIMS; j++) {
      dist2 += pow(target[j] - coords[i][j], 2);
    }
    
    if (dist2 < min_dist2) {
      min_dist2 = dist2;
      min_index = i;
    }
  }
  
  return min_index;
}


double get_dist(double * target, int index) {
  double dist2 = 0;
  for (int i = 0; i < DIP_DIMS; i++) {
    dist2 += pow(target[i] - coords[index][i], 2);
  }
  return dist2;
}


//void Cool::set_clamp_values(double * cmins, double * cmaxs) {
//  /**
//   * Set values of the clamps to use in interpolate().
//   *
//   * double * cmins        Pointer to array of doubles. Min values for clamping.
//   * double * cmaxs        Pointer to array of doubles. Max values for clamping.
//   */
//  for (int i = 0; i < D; i++) {
//    CLAMP_MAX[i] = cmaxs[i];
//    CLAMP_MIN[i] = cmins[i];
//  }
//}