//
// Created by vetinari on 14.12.20.
//

#ifndef MASTER_PROJECT_C_PART_COOLCINTERFACE_H
#define MASTER_PROJECT_C_PART_COOLCINTERFACE_H

int Cool_new();
void Cool_read_files(int, char *, char *, char *);
void Cool_reset(int);
void Cool_construct_tree(int);
void Cool_save_tree(int, char *);
double Cool_interpolate(int, double *);

int CoolManager_new(double, double, char *);
double CoolManager_interpolate(int, double *, double);
void CoolManager_save_trees(int, char *, char *);
#endif //MASTER_PROJECT_C_PART_COOLCINTERFACE_H
