//
// Created by slueders on 11/25/20.
//

#include "cool.h"
#include "CoolManager.h"
#include <map>
#include <iostream>
#include <assert.h>

#define AUTOMATIC_LOADING 1


template<int N_MAX, int D, int S_MAX>
double CoolManager<N_MAX, D, S_MAX>::interpolate(double * args, double z) {
#if AUTOMATIC_LOADING==1
    autoload(z);
#endif
    double lambda_low = low->interpolate(args);
    double lambda_high = high->interpolate(args);
    return (lambda_low * (z_high - z) + lambda_high * (z - z_low)) / z_diff;
}


template<int N_MAX, int D, int S_MAX>
void CoolManager<N_MAX, D, S_MAX>::save_trees(std::string fname_low, std::string fname_high) {
    low->save_btree(fname_low);
    high->save_btree(fname_high);
}


template<int N_MAX, int D, int S_MAX>
void CoolManager<N_MAX, D, S_MAX>::autoload(double z) {
    if (z > z_high) {
        std::cerr << "Moving backwards in time?" << std::endl;
        std::cerr << "Current z interval: [" << z_low << ", " << z_high << "]" << std::endl;
        std::cerr << "Tried to move to z = " << z << std::endl; // TODO This will lead to rpoblems if one calculation is ahead of the others
        assert(0);
    }

    if (z <= z_high && z >= z_low) {
        return;
    }

//    for (std::map<double, std::string>::iterator it=filenames.begin(); it!=filenames.end(); ++it) {
    std::string fname;
    for (auto it=filenames.begin(); it!=filenames.end(); ++it) {
        if (it->first > z_low) {
            it--;
            fname = it->second;
            break;
        }
    }

    push_slice(fname);
}


template<int N_MAX, int D, int S_MAX>
void CoolManager<N_MAX, D, S_MAX>::push_slice(std::string filename) {
    Cool<N_MAX, D, S_MAX> * temp = high;
    high = low;
    low = temp;
    low->read_files(
            filename + ".points",
            filename + ".tris",
            filename + ".neighbors"
    );
    low->construct_btree();

}

template<int N_MAX, int D, int S_MAX>
void CoolManager<N_MAX, D, S_MAX>::test() {

}