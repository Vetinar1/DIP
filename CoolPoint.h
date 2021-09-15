//
// Created by vetinari on 01.12.20.
//

#ifndef DIP_COOLPOINT_H
#define DIP_COOLPOINT_H

#include "CoolConst.h"

class Point {
    /**
     * Class representing a single data point.
     */
    friend class Cool;
    friend class Simplex;

public:
    double coords[DIP_DIMS];
    double value;
};

#endif //DIP_COOLPOINT_H
