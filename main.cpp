//
// Created by vetinari on 15.09.20.
//

#include <iostream>
#include "cool.cpp"

int main() {
    Cool<239, 2, 452> cool;

    cool.read_files("../points", "../dtriangulation", "../dneighbors");
    cool.construct_btree();

    double coord[2] = {2, 5};
    cool.interpolate(coord);

    return 0;
}