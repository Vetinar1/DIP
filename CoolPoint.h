//
// Created by vetinari on 01.12.20.
//

#ifndef MASTER_PROJECT_C_PART_COOLPOINT_H
#define MASTER_PROJECT_C_PART_COOLPOINT_H

#include "CoolConst.h"

class Point {
    /**
     * Class representing a single data point.
     */
    friend class Cool;
    friend class Simplex;
private:
    double value;

public:
    double coords[D];
};

#endif //MASTER_PROJECT_C_PART_COOLPOINT_H
