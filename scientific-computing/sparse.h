//
//  sparse.h
//  scientific-computing
//
//  Created by 马 征 on 13-11-5.
//  Copyright (c) 2013年 马 征. All rights reserved.
//

#ifndef scientific_computing_sparse_h
#define scientific_computing_sparse_h

#include "dynamicVector.h"


template<class T> struct sparseMat
//Sparse matrix data structure for compressed column storage.
{
    int nrows;      //    Number of rows.
    int ncols;      //    Number of columns.
    int nvals;      //    Maximum number of nonzeros.
    dynamicVector<int> col_ptr; //    Pointers to start of columns. Length is ncols+1.
    dynamicVector<int> row_ind; //    Row indices of nonzeros.
    dynamicVector<T> val;    //    Array of nonzero values.


    sparseMat();        //Default constructor.
    sparseMat(int m,int n,int nnvals);          //Constructor. Initializes vector to zero.
    dynamicVector<T> ax(const dynamicVector<T> &x) const;       //Multiply A by a vector x[0..ncols-1].
    dynamicVector<T> atx(const dynamicVector<T> &x) const;      //Multiply AT by a vector x[0..nrows-1].
    sparseMat transpose() const;                //Form AT .
 
};


template<class T>
sparseMat<T>::sparseMat() : nrows(0),ncols(0),nvals(0),col_ptr(), row_ind(),val() {}

template<class T>
sparseMat<T>::sparseMat(int m,int n,int nnvals) : nrows(m),ncols(n), nvals(nnvals),col_ptr(n+1,0),row_ind(nnvals,0),val(nnvals,0.0) {}

template<class T>
dynamicVector<T> sparseMat<T>::ax(const dynamicVector<T> &x) const {
    dynamicVector<T> y(nrows,0.0);
    for (int j=0;j<ncols;j++) {
        for (int i=col_ptr[j];i<col_ptr[j+1];i++)
            y(row_ind[i]) += val[i]*x[j];
    }
    return y;
}

template<class T>
dynamicVector<T> sparseMat<T>::atx(const dynamicVector<T> &x) const { dynamicVector<T> y(ncols);
    for (int i=0;i<ncols;i++) {
        y(i)=0.0;
        for (int j=col_ptr[i];j<col_ptr[i+1];j++)
            y(i) += val[j]*x[row_ind[j]];
    }
    return y;
}

template<class T>
sparseMat<T> sparseMat<T>::transpose() const {
    int i,j,k,index,m=nrows,n=ncols;
    sparseMat at(n,m,nvals);
    
    dynamicVector<int> count(m,0);
    for (i=0; i<n; i++)
        for (j=col_ptr[i];j<col_ptr[i+1];j++) {
            k=row_ind[j];
            count(k)++;
        }
    for (j=0; j<m; j++)
        at.col_ptr(j+1)=at.col_ptr[j]+count[j];
    for (j=0; j<m; j++)
        count(j)=0;
    for (i=0;i<n;i++)           //Main loop.
        for (j=col_ptr[i];j<col_ptr[i+1];j++) {
            k=row_ind[j];
            index=at.col_ptr[k]+count[k];
            at.row_ind(index)=i;
            at.val(index)=val[j];
            count(k)++;
        }
    return at;
}
























#endif
