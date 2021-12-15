//
// Created by vetinari on 01.12.20.
//

#ifndef DIP_Point_H
#define DIP_Point_H

#include "Const.h"

class Point {
    /**
     * Class representing a single data point.
     */
    friend class Cool;
    friend class Simplex;

public:
    double coords[DIP_DIMS];
    double value[DIP_VARNR];
};

#endif //DIP_Point_H
