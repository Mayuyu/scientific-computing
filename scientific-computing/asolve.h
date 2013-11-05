//
//  asolve.h
//  scientific-computing
//
//  Created by 马 征 on 13-11-5.
//  Copyright (c) 2013年 马 征. All rights reserved.
//

#ifndef scientific_computing_asolve_h
#define scientific_computing_asolve_h

#include "linbcg.h"
#include "sparse.h"
#include "dynamicVector.h"

template <class T>
struct sparseLinbcg : Linbcg<T> {
    sparseMat<T> &mat;
    int n;
    sparseLinbcg(sparseMat<T> &matrix) : mat(matrix), n(mat.nrows) {}
//    The constructor just binds a reference to your sparse matrix, making it available to asolve and atimes. To solve for a right-hand side, you call this object’s solve method, as defined in the base class.
    void atimes(dynamicVector<T> &x, dynamicVector<T> &r, const int itrnsp) {
        if (itrnsp) r=mat.atx(x);
    else r=mat.ax(x);
    }
    void asolve(dynamicVector<T> &b, dynamicVector<T> &x, const int itrnsp) {
        int i,j;
        double diag;
        for (i=0;i<n;i++) {
            diag=0.0;
            for (j=mat.col_ptr[i];j<mat.col_ptr[i+1];j++)
                if (mat.row_ind[j] == i) {
                    diag=mat.val[j];
                    break;
                }
            x(i)=(diag != 0.0 ? b[i]/diag : b[i]);
        }
    }
};

#endif
