//
// Created by vetinari on 14.12.20.
//

#ifndef DIP_COOLCINTERFACE_H
#define DIP_COOLCINTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

int Cool_new();
void Cool_read_files(int, char *, char *, char *);
void Cool_reset(int);
void Cool_construct_tree(int);
void Cool_save_tree(int, char *);
double Cool_interpolate(int, double *);
void Cool_set_clamps(int, double * mins, double * maxs);

int CoolManager_new(double, double, char *, int);
double CoolManager_interpolate(int, double *, double);
void CoolManager_save_trees(int, char *, char *);
void CoolManager_set_clamps(int, double *, double *);

#ifdef __cplusplus
};
#endif

#endif //DIP_COOLCINTERFACE_H
