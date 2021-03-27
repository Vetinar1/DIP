//
// Created by vetinari on 01.12.20.
//

#ifndef MASTER_PROJECT_C_PART_COOLSIMPLEX_H
#define MASTER_PROJECT_C_PART_COOLSIMPLEX_H

#include "CoolPoint.h"
#include "CoolConst.h"

class Simplex {
    /**
     * Class representing an D-dimensional simplex in the triangulation. Has D+1 vertices.
     *
     * int D        Number of dimensions
     */
    friend class Cool;
private:
    Point * points[D+1];             // D+1 points; Array of pointers to Point<D>
    double centroid[D];
    double midpoints[D+1][D];           // Midpoints of the faces; D+1 faces, D coordinates
    double normals[D+1][D];             // Outward pointing normals of the faces; D+1 faces, D vector components
    double T_inv[D][D];                 // T: https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Barycentric_coordinates_on_tetrahedra
//    double T[D][D];                 // T: https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Barycentric_coordinates_on_tetrahedra
    int neighbour_indices[D+1];         // One neighbour opposite every point
    Simplex * neighbour_pointers[D+1];
    double sbtree_radius_sq;            // Squared radius of ball in ball tree
    int index;

    void construct_T_inv();
    void calculate_T();
    void invert_T();
    double * convert_to_bary(const double *);
    int check_bary(const double *);
    int calculate_normals();
    void calculate_midpoints();

    template<int M, int N1, int N2> void gauss_elimination(double (&)[M][N1], double (&)[M][N2], int);
public:
    // TODO does this need to be public?
    double * find_normal(Point ** );
    void validate_simplex();
    int validate_normals();
    void print_error_info();

    // Pointers to left and right children in ball tree
    Simplex * lchild;
    Simplex * rchild;

    Simplex() {
        sbtree_radius_sq = 0;
    };

    void calculate_centroid() {
        /**
         * Calculate centroid (=avg) from points.
         */
        for (int i = 0; i < D; i++) {   // coordinates
            centroid[i] = 0;
            for (int j = 0; j < D+1; j++) {     // points
                centroid[i] += points[j]->coords[i];
            }
            centroid[i] /= (D+1);
        }
    }
};

#endif //MASTER_PROJECT_C_PART_COOLSIMPLEX_H
