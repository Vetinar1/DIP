//
// Created by Stefan Lüders on 26/07/2021.
//

#ifndef DIP_PSIBALLTREE_H
#define DIP_PSIBALLTREE_H

class PSIBallTree {
public:
    int pivot; // index in coords/vals arrays
    PSIBallTree * lchild;
    PSIBallTree * rchild;
    
    PSIBallTree * closechild;
    PSIBallTree * farchild;
    double radius;
    
    void cleanup() {
        if (lchild != nullptr) {
            lchild->cleanup();
        }
        if (rchild != nullptr) {
            rchild->cleanup();
        }
        
        delete lchild;
        delete rchild;
    };
};


#endif //DIP_PSIBALLTREE_H
