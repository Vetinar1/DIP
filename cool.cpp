//
// Created by vetinari on 15.09.20.
//

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

template<int D>
class Simplex {
    /**
     * Class representing an M-d simplex in the triangulation. Has M+1 vertices.
     *
     * int M        Number of dimensions
     */
    template<int, int, int> friend class Cool;
private:
    // TODO I would like to make these const, but I don't think I can, since the value is determined at runtime
    int points[D+1];                    // D+1 points (indices)
    int neighbour_indices[D+1];         // One neighbour oppposite every point
    Simplex * neighbour_pointers[D+1];
    double centroid[D];
    double btree_radius_sq;
public:
    Simplex * lchild;
    Simplex * rchild;

    Simplex() {
        btree_radius_sq = 0;
    };

    void calculate_centroid(const double * coords) {
        /**
         * double * coords      Pointer to array to use for lookup of actual coordinates - i.e. Cool.points
         */
        for (int i = 0; i < D; i++) {   // coordinates
            centroid[i] = 0;
            for (int j = 0; j < D+1; j++) {     // points
                centroid[i] += *(coords + points[j]*(D+1) + i);
            }
            centroid[i] /= (D+1);
        }
    }
};

template<int N, int D, int S>
class Cool {
    /**
     * int N        Number of points
     * int D        Dimensionality of points
     * int S        Number of Simplices
     */
private:
    double points[N][D+1];      // N points, D dimensions, 1 Value
    Simplex<D> simplices[S];

public:
    Cool() {};

    int read_files(std::string, std::string, std::string);
};


template<int N, int D, int S>
int Cool<N, D, S>::read_files(std::string cool_file, std::string tri_file, std::string neighbour_file) {
    /**
     * Reads files generated by the python program.
     *
     * cool_file        Path to file containing points of the grid
     * tri_file         Path to file containing triangulation
     * neighbour_file    Path to file containing triangulation neighbourhood relations
     */
    std::ifstream file;
    std::string line;
    std::string value;

    /* Read points */
    file.open(cool_file);
    if (!file.is_open()) {
        std::cerr << "Error reading " << cool_file << std::endl;
        return 1;
    }

    for (int i = 0; i < N; i++) {
        std::getline(file, line);
        std::stringstream linestream(line);
        for (int j = 0; j < D+1; j++) {     // D coordinates, 1 value
            std::getline(linestream, value, ',');
            this->points[i][j] = std::stod(value);
        }

    }

    file.close();


    /* Read Simplices */
    file.open(tri_file);
    if (!file.is_open()) {
        std::cerr << "Error reading " << tri_file << std::endl;
        return 2;
    }


    for (int i = 0; i < S; i++) {
        std::getline(file, line);
        std::stringstream linestream(line);
        for (int j = 0; j < D+1; j++) {     // D+1 points per simplex
            std::getline(linestream, value, ',');
            simplices[i].points[j] = std::stoi(value);
        }
        simplices[i].calculate_centroid(&(this->points[0][0]));
    }

    file.close();

    /* Read neighbourhood relations */
    file.open(neighbour_file);
    if (!file.is_open()) {
        std::cerr << "Error reading " << neighbour_file << std::endl;
        return 2;
    }


    for (int i = 0; i < S; i++) {
        std::getline(file, line);
        std::stringstream linestream(line);
        for (int j = 0; j < D+1; j++) {     // D+1 points per simplex
            std::getline(linestream, value, ',');
            simplices[i].neighbour_indices[j] = std::stoi(value);
            simplices[i].neighbour_pointers[j] = &(simplices[std::stoi(value)]);
        }
    }

    file.close();
    return 0;
}