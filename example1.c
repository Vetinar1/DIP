#include <stdio.h>
#include "CInterface.h"

int main() {
    printf("Interpolating using CoolManager object:\n");
    int cm_idx = CoolManager_new(3.0, 3.5, "example_data/mapfile", 0);
    double coord[3]; // T, nH, Z
    double z = 3.42;

    double clamp_mins[3] = {2.2, -8.8, -2.8};
    double clamp_maxs[3] = {8.8, 3.8, 0.8};

    CoolManager_set_clamps(cm_idx, clamp_mins, clamp_maxs);

    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            for (int k = 0; k < 10; k++) {
                coord[0] = clamp_mins[0] + i * (clamp_maxs[0]-clamp_mins[0])/10.;
                coord[1] = clamp_mins[1] + j * (clamp_maxs[1]-clamp_mins[1])/10.;
                coord[2] = clamp_mins[2] + k * (clamp_maxs[2]-clamp_mins[2])/10.;
                CoolManager_interpolate(cm_idx, &(coord[0]), z);
            }
        }
    }


    printf("\n\n\n\n\n\nInterpolating using Cool object:\n\n\n\n");
    int c_idx = Cool_new(0);
    Cool_read_files(c_idx, "example_data/z3.0.points", "example_data/z3.0.tris", "example_data/z3.0.neighbors");
    Cool_construct_tree(c_idx);
    Cool_set_clamps(c_idx, clamp_mins, clamp_maxs);

    for (int i = 0; i < 100; i++) {
        for (int j = 0; j < 100; j++) {
            coord[0] = clamp_mins[0] + i * (clamp_maxs[0]-clamp_mins[0])/100.;
            coord[1] = clamp_mins[1] + j * (clamp_maxs[1]-clamp_mins[1])/100.;
            double interp = Cool_interpolate(c_idx, &(coord[0]));
        }
    }
    return 0;
}
