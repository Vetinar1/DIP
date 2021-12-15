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

    PSIBallTree() {
	lchild = nullptr;
	rchild = nullptr;

	closechild = nullptr;
	farchild = nullptr;
    }
    
    void cleanup() {
        if (lchild != nullptr) {
            lchild->cleanup();
	    delete lchild;
        }
        if (rchild != nullptr) {
            rchild->cleanup();
	    delete rchild;
        }
    };
};


#endif //DIP_PSIBALLTREE_H
