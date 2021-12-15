//
// Created by vetinari on 01.12.20.
//

#ifndef DIP_Simplex_H
#define DIP_Simplex_H

#include "Point.h"
#include "Const.h"

class Simplex {
    /**
     * Class representing an DIP_DIM-dimensional simplex in the triangulation. Has DIP_DIM+1 vertices.
     *
     * int DIP_DIM        Number of dimensions
     */
    friend class Cool;
private:
    double centroid[DIP_DIMS];
    double midpoints[DIP_DIMS+1][DIP_DIMS];           // Midpoints of the faces; DIP_DIM+1 faces, DIP_DIM coordinates
    double normals[DIP_DIMS+1][DIP_DIMS];             // Outward pointing normals of the faces; DIP_DIM+1 faces, DIP_DIM vector components
    double T_inv[DIP_DIMS][DIP_DIMS];                 // T: https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Barycentric_coordinates_on_tetrahedra
//    double T[DIP_DIMS][DIP_DIMS];                 // T: https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Barycentric_coordinates_on_tetrahedra
    int neighbor_indices[DIP_DIMS+1];         // One neighbor opposite every point
    Simplex * neighbor_pointers[DIP_DIMS+1];
    double sbtree_radius_sq;            // Squared radius of ball in ball tree
    int index;
    
    void calculate_T();
    void invert_T();
    int calculate_normals();
    void calculate_midpoints();

    template<int M, int N1, int N2> void gauss_elimination(double (&)[M][N1], double (&)[M][N2], int);
public:
    Point * points[DIP_DIMS+1];             // DIP_DIM+1 points; Array of pointers to Point<DIP_DIM>
    void construct_T_inv();
    double * convert_to_bary(const double *);
    int check_bary(const double *);
    
    double * find_normal(Point ** );
    void validate_simplex();
    int validate_normals();
    void print_error_info();
    
    double get_quality();

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
        for (int i = 0; i < DIP_DIMS; i++) {   // coordinates
            centroid[i] = 0;
            for (int j = 0; j < DIP_DIMS+1; j++) {     // points
                centroid[i] += points[j]->coords[i];
            }
            centroid[i] /= (DIP_DIMS+1);
        }
    }
};

#endif //DIP_Simplex_H
