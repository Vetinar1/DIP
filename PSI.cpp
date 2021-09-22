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


double get_coord(int i, int j) {
  return coords[i][j];
}


double get_val(int i, int j) {
  return vals[i][j];
}


double dot(double * v, double * w) {
  /**
   * Simple dot product between two vectors.
   */
  double out = 0;
  for (int i = 0; i < DIP_DIMS; i++) {
    out += v[i] * w[i];
  }
  return out;
}



double get_dist(double * target, int index) {
  /**
   * Helper function. Get distance from target to point with given index.
   */
  double dist2 = 0;
  for (int i = 0; i < DIP_DIMS; i++) {
    dist2 += pow(target[i] - coords[index][i], 2);
  }
  return dist2;
}


void psi_init() {
  /**
   * Set up global arrays for interpolation
   */
  coords = new double * [DIP_NMAX];
  vals   = new double * [DIP_NMAX];
  
  for (int i = 0; i < DIP_NMAX; i++) {
    coords[i] = new double[DIP_DIMS];
    vals[i]   = new double[DIP_VARNR];
  }
}


int psi_read_points(std::string cool_file) {
  /**
   * Read data from file, see also equivalent in CoolCool.cpp
   *
   * @param cool_file   .points file returned from CHIPS
   */
  std::ifstream file;
  std::string line;
  std::string value;
  
  
  /* Read points */
  file.open(cool_file);
  if (!file.is_open()) {
    std::cerr << "Error reading " << cool_file << std::endl;
    return 1;
  } else {
    std::cout << "Reading " << cool_file << " (max " << DIP_NMAX << " lines)" << std::endl;
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
  
  psi_construct_btree(coords);
  
  return 0;
}


PSIBallTree * psi_construct_btree_recursive(double ** base, int * indices, int n) {
  /**
   * Recursively construct btree
   *
   * Do not call this function directly. Use construct_btree() instead.
   *
   * @param base        The array containing *all* points that the tree is being built on. Actual coordinate values
   * @param indices     The array containing the points that the *current subtree* is being built on. Indices of base
   * @param n           Size of array indices
   * @return            Ball tree node that contains all the points in indices
   */
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
      coord_i[j] = base[indices[j]][i];
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
    double val = base[indices[i]][lspread_dim];
    double dist = fabs(median - val);
    
    if (dist < min_dist) {
      // Add old pivot to L or R
      if (node->pivot != -1) {
        if (base[node->pivot][lspread_dim] <= median) {
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
    node->lchild = psi_construct_btree_recursive(base, L, lcount);
  } else {
    node->lchild = nullptr;
  }
  if (rcount > 0) {
    node->rchild = psi_construct_btree_recursive(base, R, rcount);
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
      dist += pow(base[indices[i]][j] - base[node->pivot][j], 2);
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
      lchild_dist2 += pow(base[node->pivot][i] - base[node->lchild->pivot][i], 2);
      rchild_dist2 += pow(base[node->pivot][i] - base[node->rchild->pivot][i], 2);
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


int psi_construct_btree(double ** points) {
  /**
   * Builds a ball tree for quick nearest neighbour finding. Public adapter for
   * construct_btree_recursive.
   *
   * @param points      Array of points to build balltree on
   */
  // Construct array of indices
  int * indices = new int[N_LIM];
  for (int i = 0; i < N_LIM; i++) {
    indices[i] = i;
  }
  
  btree = psi_construct_btree_recursive(points, indices, N_LIM);
  delete[] indices;
  
  return 0;
}


void psi_find_k_nearest_neighbour_recursive(PSIBallTree * root, const double * target, double ** base,
                                            std::priority_queue<distpoint> * Q, int k) {
  /**
   * Efficiently find the k nearest neighbours of target in coords.
   * This function may look a bit odd because it was originally only supposed to find *one* nearest neighbour.
   * Do not use this function directly, use psi_find_k_nearest_neighbour() instead.
   *
   * @param root        Ball tree to search in
   * @param target      Pointer to coordinates of target
   * @param base        Array of coordinates that the ball tree was built on
   * @param Q           Priority queue containing the currently known k nearest neighbours
   * @param k           How many nearest neighbours to find
   *
   * Returns nothing, the output is the priority queue, which is modified in-place.
   */
  double root_dist = 0;
  double Q_dist = 0;

  for (int i = 0; i < DIP_DIMS; i++) {
    root_dist += pow(target[i] - base[root->pivot][i], 2);
  }
  root_dist = sqrt(root_dist);
  
  if (!Q->empty()) {
    // Distance to current node and distance to farthest node in Q
    distpoint dp = Q->top();
  
    for (int i = 0; i < DIP_DIMS; i++) {
      Q_dist += pow(target[i] - base[dp.index][i], 2);
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
    psi_find_k_nearest_neighbour_recursive(root->closechild, target, base, Q, k);
  }
  if (root->farchild != nullptr) {
    psi_find_k_nearest_neighbour_recursive(root->farchild, target, base, Q, k);
  }
}


int * psi_find_k_nearest_neighbor(double * target, double ** points, PSIBallTree * btree, int k) {
  /**
   * Public adapter for psi_find_k_nearest_neighbour_recursive(). Efficently finds the k nearest neighbours of
   * target in points.
   *
   * @param target  Pointer to the target coordinates
   * @param points  Array of points to find nearest neighbour in
   * @param btree   Balltree built on points
   * @param k       How many nearest neighbours to find
   * @return        Pointer to array of indices of k nearest neighbours in coords
   */
  std::priority_queue<distpoint> * Q = new std::priority_queue<distpoint>;
  psi_find_k_nearest_neighbour_recursive(btree, target, points, Q, k);
  
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


int * psi_find_k_nearest_neighbor(double * target, int k) {
  /**
   * Find k nearest neighbours in coords
   */
  return psi_find_k_nearest_neighbor(target, coords, btree, k);
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
    // TODO try using kd tree after all...
    int nn;
    if (d == DIP_DIMS) {
      nn = 0;
    } else {
#if 0
      int * indices = new int[k];
      int n = 0;
      for (int i = 0; i < k; i++) {
        if (neigh_mask[i]) {
          indices[n] = i;
          n++;
        }
      }
      double min_dist2 = DBL_MAX;
      PSIBallTree * best = nullptr;
      PSIBallTree * temptree = psi_construct_btree_recursive(neigh_coords, indices, n);
      PSIBallTree * result = psi_find_nearest_neighbour_recursive(temptree, ptarget, neigh_coords, best, &min_dist2);
      nn = result->pivot;
#else
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
#endif
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
    if (n_count < 1 + d) {
#ifdef PSI_SHOW_DIAGNOSTICS
      std::cout << "n_count is " << n_count << " at d = " << d << std::endl;
      std::cout << std::endl;
#endif
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
    double linevec[DIP_DIMS];
    double * proj1d = new double[k];
    int linevecflag = 0;
    
    double posmin =  DBL_MAX;
    double negmax = -DBL_MAX;
    int posind = -1;
    int negind = -1;
  
    for (int i = 0; i < k; i++) {
      if (neigh_mask[i] == 0) {
        continue;
      }
    
      if (!linevecflag) {
        // Calculate vector parallel to line
        for (int j = 0; j < DIP_DIMS; j++) {
          linevec[j] = neigh_coords[i][j] - ptarget[j];
        }
      
        linevecflag = 1;
      }
    
      // Project difference vectors between other points and ptarget on linevec
      double tempvec[DIP_DIMS];
      for (int j = 0; j < DIP_DIMS; j++) {
        tempvec[j] = neigh_coords[i][j] - ptarget[j];
      }
    
      proj1d[i] = dot(linevec, tempvec);
  
      if (proj1d[i] > 0 && proj1d[i] < posmin) {
        posmin = proj1d[i];
        posind = i;
      } else if (proj1d[i] < 0 && proj1d[i] > negmax) {
        negmax = proj1d[i];
        negind = i;
      }
    }
  
    // If all scalar products were positive or negative, ptarget is an "endpoint" on the line and there is no solution
    if (posind == -1 || negind == -1) {
#ifdef PSI_SHOW_ERRORS
      std::cout << "Couldn't find two points on line; posind " << posind << " - negind " << negind << std::endl;
      int sum = 0;
      for (int i = 0; i < k; i++) {
        sum += neigh_mask[i];
      }
      std::cout << sum << " points were left" << std::endl;
#endif
      delete[] simplex;
      simplex = nullptr;
    } else {
      simplex[1] = neighbours[negind];
      simplex[0] = neighbours[posind];
    }
  
    delete[] proj1d;
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
  adaptive_calls++;
  adaptive_executions++;
  int * neighbours = psi_find_k_nearest_neighbor(target, k);
  int * simplex = psi_projective_simplex_algorithm(neighbours, target, k);
  
  int iterations = 0;
  while (simplex == nullptr && iterations < max_steps) {
    k = (int) k * factor;
    neighbours = psi_find_k_nearest_neighbor(target, k);
    simplex = psi_projective_simplex_algorithm(neighbours, target, k);
    iterations++;
    adaptive_executions++;
  }

#ifdef PSI_SHOW_DIAGNOSTICS
  if (iterations) {
    std::cout << "Extra iterations: " << iterations << std::endl;
  }
#endif
  
  return simplex;
}

