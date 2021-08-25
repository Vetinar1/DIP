//
// Created by Stefan Lüders on 26/07/2021.
//

#include "CoolSimplex.h"
#include "CoolPoint.h"
#include "CoolConst.h"
#include "PSIBallTree.h"
#include "PSI.h"
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


int get_nlim() {
  return N_LIM;
}

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
    leaf->closechild = nullptr;
    leaf->farchild = nullptr;
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
    assert(node->lchild != nullptr || node->rchild != nullptr);
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


int psi_find_nearest_neighbour(double * target) {
  double min_dist2 = DBL_MAX;
  PSIBallTree * best = nullptr;
  PSIBallTree * node = psi_find_nearest_neighbour_recursive(btree, target, best, &min_dist2);
  
  return node->pivot;
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


void psi_find_k_nearest_neighbour_recursive(PSIBallTree * root, const double * target,
                                            std::priority_queue<distpoint> * Q, int k) {
  /**
   * This function may look a bit odd because it was originally only supposed to find *one* nearest neighbour.
   */
  double root_dist = 0;
  double Q_dist = 0;

  for (int i = 0; i < DIP_DIMS; i++) {
    root_dist += pow(target[i] - coords[root->pivot][i], 2);
  }
  root_dist = sqrt(root_dist);
  
  if (!Q->empty()) {
    // Distance to current node and distance to farthest node in Q
    distpoint dp = Q->top();
  
    for (int i = 0; i < DIP_DIMS; i++) {
      Q_dist += pow(target[i] - coords[dp.index][i], 2);
    }
    Q_dist = sqrt(Q_dist);
  
    // Recursion exit condition - target point is further outside the current node's ball than the distance
    // to the current farthest nearest neighbor -> there can't be a closer neighbor in this ball/subtree
    if (root_dist - root->radius >= Q_dist) {
      return;
    }
  }
  
  // The current node is closer than at least one node from the queue. Add it, keep queue length <= k
  Q->push(distpoint{root->pivot, root_dist});
  if (Q->size() > k) {
    Q->pop();
  }
  
  if (root->closechild != nullptr) {
    psi_find_k_nearest_neighbour_recursive(root->closechild, target, Q, k);
  }
  if (root->farchild != nullptr) {
    psi_find_k_nearest_neighbour_recursive(root->farchild, target, Q, k);
  }
}


int * psi_find_k_nearest_neighbor(double * target, int k) {
  std::priority_queue<distpoint> * Q = new std::priority_queue<distpoint>;
  psi_find_k_nearest_neighbour_recursive(btree, target, Q, k);
  
  assert(Q->size() == k);
  
  int * out = new int[k];
  for (int i = k-1; i >= 0; i--) {
    distpoint temp = Q->top();
    out[i] = temp.index;
    Q->pop();
  }
  
  assert(Q->size() == 0);
  delete Q;

  return out;
}


int * psi_find_k_nearest_neighbour_bruteforce(double * target, int k) {
  std::priority_queue<distpoint> Q;
  for (int i = 0; i < N_LIM; i++) {
    double dist2 = get_dist(target, i);
    if (dist2 < DIP_EPSILON) {
      continue;
    }
    Q.push(distpoint{i, dist2});
    
    if (Q.size() > k) {
      Q.pop();
    }
  }
  
  int * out = new int[k];
  for (int i = k-1; i >= 0; i--) {
    distpoint temp = Q.top();
    out[i] = temp.index;
    Q.pop();
  }
  
  assert(Q.empty());
  
  return out;
}


double get_dist(double * target, int index) {
  double dist2 = 0;
  for (int i = 0; i < DIP_DIMS; i++) {
    dist2 += pow(target[i] - coords[index][i], 2);
  }
  return dist2;
}

double dot(double * v, double * w) {
  double out = 0;
  for (int i = 0; i < DIP_DIMS; i++) {
    out += v[i] * w[i];
  }
  return out;
}


int * psi_projective_simplex_algorithm(int * neighbours, double * target, int k) {
  // Copy relevant neighbours into separate array of dimensions k x DIP_DIMS
  // The DIP_DIMS columns are going to hold the actual coordinates; these are going to be
  // modified as we do the projection in each iteration
  double ** neigh_coords = new double * [k];
  for (int i = 0; i < k; i++) {
    neigh_coords[i] = new double[DIP_DIMS];
    for (int j = 0; j < DIP_DIMS; j++) {
      neigh_coords[i][j] = coords[neighbours[i]][j];
    }
  }
  
  // As we filter out neighbours during the projection process we will not want to iterate over the entire array
  // Use this array as a mask to keep track of which neighbours are still in play
  int * neigh_mask = new int[k];
  for (int i = 0; i < k; i++) {
    neigh_mask[i] = 1;
  }
  
  int * simplex = new int[DIP_DIMS+1];  // The indices of the points that are going to form the simplex
  double diff[DIP_DIMS];                // difference vector target -> nn; defines projection plane
  double ptarget[DIP_DIMS];             // Projected target vector
  
  for (int i = 0; i < DIP_DIMS; i++) {
    ptarget[i] = target[i];
  }
  
  
  int d = DIP_DIMS;
  int failflag = 0;
  while (d > 1) {
    // Find nearest neighbor. In the first iteration this is always the first neighbour in the list, otherwise
    // do a brute force search
    int nn;
    if (d == DIP_DIMS) {
      nn = 0;
    } else {
      double min_dist2 = DBL_MAX;
      for (int i = 0; i < k; i++) {
        if (neigh_mask[i] == 0) {
          continue;
        }
        
        double dist2 = 0;
        for (int j = 0; j < DIP_DIMS; j++) {
          dist2 += pow(ptarget[j] - neigh_coords[i][j], 2);
        }
        
        if (dist2 < min_dist2) {
          min_dist2 = dist2;
          nn = i;
        }
      }
    }
    
    // Add the nearest neighbor to the solution
    simplex[d] = neighbours[nn];
    // Remove it from further considerations
    neigh_mask[nn] = 0;
    
    // calculate difference vector (and its squared length);
    double difflen2 = 0;
    for (int i = 0; i < DIP_DIMS; i++) {
      diff[i] = neigh_coords[nn][i] - ptarget[i];
      difflen2 += diff[i] * diff[i];
    }
    
    
    // project and filter
    for (int i = 0; i < k; i++) {
      if (neigh_mask[i] == 0) {
        continue;
      }
      
      // 1. project on normal
      double pn[DIP_DIMS];        // projection on normal
      double pn_shift[DIP_DIMS];  // projection on normal, shifted to ptarget
      double pn_factor = dot(neigh_coords[i], diff) / difflen2;
      for (int j = 0; j < DIP_DIMS; j++) {
        pn[j]       = pn_factor * diff[j];
        pn_shift[j] = pn[j] - ptarget[j];
//        std::cout << pn[j] << "\t" << pn_shift[j] << "\t" << diff[j] << std::endl;
      }
      
      // 2. Filter (only keep neighbours on "negative" side of the plane)
      if (dot(pn_shift, diff) > 0) {
//        std::cout << "d = " << d << ": setting " << i << " to 0" << std::endl;
//        std::cout << (dot(pn_shift, diff) > 0) << std::endl;
        neigh_mask[i] = 0;
        continue;
      }
      
      // 3. Finish projection on to plane
      for (int j = 0; j < DIP_DIMS; j++) {
        neigh_coords[i][j] = neigh_coords[i][j] - pn[j];
      }
    }
    
    // Verify we still have enough neighbours to continue next iteration. We need at least 2 at the end (on the line)
    int n_count = 0;
    for (int i = 0; i < k; i++) {
      n_count += neigh_mask[i];
    }
    if (n_count < 2) {
      std::cout << "n_count is " << n_count << " at d = " << d << std::endl;
      std::cout << std::endl;
      failflag = 1;
      delete[] simplex;
      simplex = nullptr;
      break;
    }
    
    // Update ptarget (= project ptarget onto plane)
    double pt_factor = dot(ptarget, diff) / difflen2;
    for (int i = 0; i < DIP_DIMS; i++) {
      ptarget[i] = ptarget[i] - pt_factor * diff[i];
    }
    
    d--;
  }
  
  if (!failflag) {
    // Now that the iteration is complete, all points are projected onto a line. Find the closest point in both
    // directions on that line.
    // 1. Find nearest neighbor
    int nn;
    double min_dist2 = DBL_MAX;
    for (int i = 0; i < k; i++) {
      if (neigh_mask[i] == 0) {
        continue;
      }
      
      for (int j = 0; j < DIP_DIMS; j++) {
        std::cout << neigh_coords[i][j] << " ";
      }
      std::cout << std::endl;
    
      double dist2 = 0;
      for (int j = 0; j < DIP_DIMS; j++) {
        dist2 += pow(ptarget[i] - neigh_coords[i][j], 2);
      }
    
      if (dist2 < min_dist2) {
        min_dist2 = dist2;
        nn = i;
      }
    }
    std::cout << "ptarget: " << ptarget[0] << " " << ptarget[1] << " " << ptarget[2] << std::endl;
  
    // simplex was filled backwards, so the second to last element to fill is 1
    simplex[1] = neighbours[nn];
    neigh_mask[nn] = 0;
    
//    std::cout << "simplex so far" << std::endl;
//    for (int i = DIP_DIMS; i >= 0; i--) {
//      std::cout << simplex[i] << std::endl;
//    }
    
    // 2. Difference vector
    double diff[DIP_DIMS];
    for (int i = 0; i < DIP_DIMS; i++) {
      diff[i] = neigh_coords[nn][i] - ptarget[i];
    }
    
    // 3. Find nearest neighbor again, but only consider those pointing in the opposite direction of diff
    min_dist2 = DBL_MAX;
    for (int i = 0; i < k; i++) {
      if (neigh_mask[i] == 0) {
        continue;
      }
//      std::cout << "i: " << i << std::endl;
      
      double finaldiff[DIP_DIMS];
      for (int j = 0; j < DIP_DIMS; j++) {
        finaldiff[j] = neigh_coords[i][j] - ptarget[j];
//        std::cout << diff[j] << " " << finaldiff[j] << std::endl;
      }
      
      if (dot(diff, finaldiff) > 0) { // same direction
//        std::cout << "wrong direction" << std::endl;
        continue;
      }
    
      double dist2 = 0;
      for (int j = 0; j < DIP_DIMS; j++) {
        dist2 += pow(finaldiff[j], 2);
      }
    
      if (dist2 < min_dist2) {
        min_dist2 = dist2;
        nn = i;
      }
    }
    
    if (min_dist2 == DBL_MAX) {
      std::cout << "couldnt find two points on line out of ";
      int sum = 0;
      for (int i = 0; i < k; i++) {
        sum += neigh_mask[i];
      }
      std::cout << sum << std::endl;
      delete[] simplex;
      simplex = nullptr;
    } else {
      simplex[0] = neighbours[nn];
    }
  }
  
  
  // Cleanup
  for (int i = 0; i < k; i++) {
    delete[] neigh_coords[i];
  }
  delete[] neigh_coords;
  delete[] neigh_mask;
  
  
  return simplex;
}


int * psi_adaptive_projective_simplex_algorithm(double * target, int k, double factor, int max_steps) {
  int * neighbours = psi_find_k_nearest_neighbor(target, k);
  int * simplex = psi_projective_simplex_algorithm(neighbours, target, k);
  
  int iterations = 0;
  while (simplex == nullptr && iterations < max_steps) {
    k = (int) k * factor;
    neighbours = psi_find_k_nearest_neighbor(target, k);
    simplex = psi_projective_simplex_algorithm(neighbours, target, k);
    iterations++;
  }
  
  if (iterations) {
    std::cout << "Extra iterations: " << iterations << std::endl;
  }

  return simplex;
}


double get_coord(int i, int j) {
  return coords[i][j];
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