#include <stdio.h>
#include "CoolCInterface.h"

#define CM
int main() {
#ifdef CM
    int cm_idx = CoolManager_new(7.5, 8, "cooldata3/mapfile", 0);
    double coord[4]; // T, nH
    double z = 7.76;
//    double clamps_min[4] = {2, -6, -5, 6};
    double clamps_min[4] = {2.2, -5.8, -4.8, 6.2};
//    double clamps_max[4] = {9, 2, 3, 12};
    double clamps_max[4] = {8.8, 1.8, 2.8, 11.8};
    CoolManager_set_clamps(cm_idx, &(clamps_min[0]), &(clamps_max[0]));
    int cnt = 0;
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            for (int k = 0; k < 10; k++) {
                for (int l = 0; l < 10; l++) {

                    coord[0] = clamps_min[0] + i * (clamps_max[0]-clamps_min[0])/10.;
                    coord[1] = clamps_min[1] + j * (clamps_max[1]-clamps_min[1])/10.;
                    coord[2] = clamps_min[2] + k * (clamps_max[2]-clamps_min[2])/10.;
                    coord[3] = clamps_min[3] + l * (clamps_max[3]-clamps_min[3])/10.;
                    printf("cnt %d\n", cnt);
                    cnt++;
                    CoolManager_interpolate(cm_idx, &(coord[0]), z);

                }
            }
        }
    }
//    double coord[2];
//    coord[1] = 3;
//    double z = 4;
//    for (int i = 0; i < 100; i++) {
//        for (int j = 0; j < 100; j++) {
//            printf("%i %f\n", i, z);
//            coord[0] = 2 + j * (8-2)/100.;
////            z = 4 - i * 4 / 100.;
//            z = 8.8;
//            double interp = CoolManager_interpolate(cm_idx, &(coord[0]), z);
//        }
//    }
    return 0;
#endif


#ifdef C
    int c_idx = Cool_new();
    Cool_read_files(c_idx, "data2d/data.csv", "data2d/dtri.csv", "data2d/dneighbours.csv");
//    printf("Reset\n");
//    Cool_reset();
    printf("Tree\n");
    Cool_construct_tree(c_idx);

    printf("Interpolate\n");
    double coord[2];
    for (int i = 0; i < 100; i++) {
        for (int j = 0; j < 100; j++) {
            printf("%i %i\n", i, j);
            coord[0] = 2 + i * (8-2)/100.;
            coord[1] = -4 + j * 8 / 100.;
            double interp = Cool_interpolate(c_idx, &(coord[0]));
        }
    }
    return 0;
#endif
}
