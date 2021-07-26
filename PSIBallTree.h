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
    double radius_sq;
    
//    void free();
};


//void PSIBallTree::free() {
//  if (lchild != nullptr) {
//    lchild->free();
//    delete lchild;
//  }
//
//  if (rchild != nullptr) {
//    rchild->free();
//    delete rchild;
//  }
//}
#endif //DIP_PSIBALLTREE_H
