#include <stdio.h>
#include "CoolCInterface.h"

#define CM
int main() {
#ifdef CM
    CoolManager_new(3.9, 4.0, "slice3d/mapfile");
    double coord[2];
    coord[1] = 3;
    double z = 4;
    for (int i = 0; i < 100; i++) {
        for (int j = 0; j < 100; j++) {
            printf("%i %f\n", i, z);
            coord[0] = 2 + j * (8-2)/100.;
            z = 4 - i * 4 / 100.;
            double interp = CoolManager_interpolate(&(coord[0]), z);
        }
    }
    return 0;
#endif


#ifdef C
    Cool_new();
    Cool_read_files("data2d/data.csv", "data2d/dtri.csv", "data2d/dneighbours.csv");
//    printf("Reset\n");
//    Cool_reset();
    printf("Tree\n");
    Cool_construct_tree();

    printf("Interpolate\n");
    double coord[2];
    for (int i = 0; i < 100; i++) {
        for (int j = 0; j < 100; j++) {
            printf("%i %i\n", i, j);
            coord[0] = 2 + i * (8-2)/100.;
            coord[1] = -4 + j * 8 / 100.;
            double interp = Cool_interpolate(&(coord[0]));
        }
    }
    return 0;
#endif
}
