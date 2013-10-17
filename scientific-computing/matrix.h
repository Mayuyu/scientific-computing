//
//  matrix.h
//  scientific-computing
//
//  Created by 马 征 on 13-10-17.
//  Copyright (c) 2013年 马 征. All rights reserved.
//

#ifndef scientific_computing_matrix_h
#define scientific_computing_matrix_h

#include "vector.h"

typedef vector<double,2> point;

typedef vector<double,3> point3;

template<class T, int N, int M> class matrix : public vector<vector<T,N>,M>{
public:
    matrix(){}
    matrix(const vector<vector<T,N>,M>&){}
    matrix(double d){
        this->set(0,point(d,0.));
        this->set(1,point(0.,d));
    }  //  constructor
    matrix(const vector<T,N>&u, const vector<T,N>&v){
        set(0,u);
        set(1,v);
    }  //  constructor
    matrix(const vector<T,N>&u, const vector<T,N>&v, const vector<T,N>&w){
        set(0,u);
        set(1,v);
        set(2,w);
    }  //  constructor
    ~matrix(){}  //  destructor
    vector<T,N>& operator()(int i){
        return vector<vector<T,N>,M>::operator()(i);
    }  //  read/write ith column
    const T& operator()(int i,int j) const{return (*this)[j][i];}//A(i,j)
    T& operator()(int i,int j, const char*){
        return (*this)(j)(i);
    }  //  read/write A(i,j)
    const matrix& operator+=(const matrix&);
    const matrix& operator-=(const matrix&);
    const matrix& operator*=(const T&);
    const matrix& operator/=(const T&);
};

template<class T, int N, int M>
ostream& operator <<(ostream& out, const matrix<T, N, M>&v){
    //    for(int i = 0;i < M; i++){
    //        for(int j = 0;j < N; j++)
    //            printf("v[%d,%d]=%f;  ",i,j,v(i,j));
    //        printf("\n");
    //    }
    
}  //  printing a matrix

template<class T, int N, int M>
const matrix<T,N,M>& matrix<T,N,M>::operator+=(const matrix<T,N,M>&m){
    vector<vector<T,N>,M>::operator+=(m);
    return *this;
}  //  adding a matrix

template<class T, int N, int M>
const matrix<T,N,M>& matrix<T,N,M>::operator-=(const matrix<T,N,M>&m){
    vector<vector<T,N>,M>::operator-=(m);
    return *this;
}  //  subtracting a matrix

template<class T, int N, int M>
const matrix<T,N,M>& matrix<T,N,M>::operator*=(const T&a){
    for(int i=0; i<M; i++)
        set(i,(*this)[i] * a);
    return *this;
}  //  multiplication by scalar

template<class T, int N, int M>
const matrix<T,N,M>& matrix<T,N,M>::operator/=(const T&a){
    for(int i=0; i<M; i++)
        set(i,(*this)[i] / a);
    return *this;
}  //  division by scalar

template<class T, int N, int M>
const matrix<T,N,M> operator*(const T&a,const matrix<T,N,M>&m){
    return matrix<T,N,M>(m) *= a;
}  //  scalar times matrix

template<class T, int N, int M>
const matrix<T,N,M> operator*(const matrix<T,N,M>&m, const T&a){
    return matrix<T,N,M>(m) *= a;
}  //  matrix times scalar

template<class T, int N, int M>
const matrix<T,N,M> operator/(const matrix<T,N,M>&m, const T&a){
    return matrix<T,N,M>(m) /= a;
}  //  matrix divided by scalar

template<class T, int N, int M>
const matrix<T,N,M> operator+(const matrix<T,N,M>&m1, const matrix<T,N,M>&m2){
    return matrix<T,N,M>(m1) += m2;
}  //  matrix plus matrix

template<class T, int N, int M>
const matrix<T,N,M> operator-(const matrix<T,N,M>&m1, const matrix<T,N,M>&m2){
    return matrix<T,N,M>(m1) -= m2;
}  //  matrix minus matrix

template<class T, int N, int M> const vector<T,M>
operator*(const vector<T,N>&v,const matrix<T,N,M>&m){
    vector<T,M> result;
    for(int i=0; i<M; i++)
        result.set(i, v * m[i]);
    return result;
}  //  vector times matrix

template<class T, int N, int M> const vector<T,N>
operator*(const matrix<T,N,M>&m,const vector<T,M>&v){
    vector<T,N> result;
    for(int i=0; i<M; i++)
        result += v[i] * m[i];
    return result;
}  //  matrix times vector

template<class T, int N, int M, int K> const matrix<T,N,K>
operator*(const matrix<T,N,M>&m1,const matrix<T,M,K>&m2){
    matrix<T,N,K> result;
    for(int i=0; i<K; i++)
        result.set(i,m1 * m2[i]);
    return result;
}  //  matrix times matrix

template<class T, int N, int M, int K> const matrix<T,N,K>&
operator*=(matrix<T,N,M>&m1,const matrix<T,M,K>&m2){
    return m1 = m1 * m2;
}  //  matrix times matrix

template<class T, int N, int M>
const matrix<T, M, N>
transpose(const matrix<T, N, M>&A){
    matrix<T, M, N> At;
    for(int i=0; i<N; i++)
        for(int j=0; j<M; j++)
            At(j,i,"write") = A(i,j);
    return At;
}  //  transpose of a matrix

template<class T, int N>
void diagonalize(matrix<T,N,N>&W,
                 const matrix<T,N,N>&A,
                 matrix<T,N,N>&Z, vector<T,N>&D){
    for(int i=0; i<N; i++){
        W(i,i,"write") = 1.;
        Z(i,i,"write") = 1.;
    }
    matrix<T,N,N> At = transpose(A);
    for(int i=0; i<N; i++){
        vector<T,N> r = At[i];
        for(int j=0; j<i; j++)
            W(i) -= r * Z[j] / D[j] * W[j];
        vector<T,N> c = A[i];
        for(int j=0; j<i; j++)
            Z(i) -= W[j] * c / D[j] * Z[j];
        D(i) = W[i] * (A * Z[i]);
    }
}  //  diagonalize a matrix

template<class T, int N>
const T det(const vector<T,N>&D){
    T result = 1.;
    for(int i=0; i<N; i++)
        result *= D[i];
    return result;
}  //  determinant of a diagonalized matrix


#endif
