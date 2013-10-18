//
//  dynamicVector.h
//  scientific-computing
//
//  Created by 马 征 on 13-10-17.
//  Copyright (c) 2013年 马 征. All rights reserved.
//

#ifndef scientific_computing_dynamicVector_h
#define scientific_computing_dynamicVector_h

#include <iostream>
#include <cmath>
#include "complex.h"

using namespace std;

int power(int basis, unsigned exp) {
    return exp ? basis * power(basis,exp-1) : 1;
}  //  "basis" to the "exp"

template <class T>
const T fabs(const T&a) {
    return a > 0 ? a: -a;
}



template<class T> class dynamicVector{
protected:
    int dimension;
    T* component;
public:
    dynamicVector(int = 0, const T& = 0);
    dynamicVector(const dynamicVector&);
    dynamicVector(const dynamicVector&, const dynamicVector&);
    const dynamicVector& operator=(const dynamicVector&);
    const dynamicVector& operator=(const T&);
    ~dynamicVector(){delete [] component;}  //  destructor
    int dim() const{ return dimension; }  //  return the dimension
    T& operator()(int i){ return component[i]; }  //read/write ith component
    const T& operator[](int i) const{ return component[i]; }  //read only
    const dynamicVector& operator+=(const dynamicVector&);
    const dynamicVector& operator-=(const dynamicVector&);
    const dynamicVector& operator*=(const T&);
    const dynamicVector& operator/=(const T&);
};

template<class T>
dynamicVector<T>::dynamicVector(int dim, const T& a) : dimension(dim),
component(dim ? new T[dim] : 0){
    for(int i = 0; i < dim; i++)
        component[i] = a;
}  //  constructor

template<class T>
dynamicVector<T>::dynamicVector(const dynamicVector<T>& v) : dimension(v.dimension),
component(v.dimension ? new T[v.dimension] : 0){
    for(int i = 0; i < v.dimension; i++)
        component[i] = v.component[i];
}  //  copy constructor

template<class T>
dynamicVector<T>::dynamicVector(const dynamicVector<T>& u,const dynamicVector<T>& v)
: dimension(u.dimension+v.dimension), component((u.dimension+v.dimension) ? new T[u.dimension+v.dimension] : 0){
    for(int i = 0; i < u.dimension; i++)
        component[i] = u.component[i];
    for(int i = 0; i < v.dimension; i++)
        component[i+u.dimension] = v.component[i];
    //    Cons2++;
    //    Dims2 += dimension;
}  //  constructor from 2 vectors

template<class T>
const dynamicVector<T>& dynamicVector<T>::operator=(const dynamicVector<T>& v){
    if(this != &v){
        if(dimension > v.dimension)
            delete [] (component + v.dimension);
        if(dimension < v.dimension){
            delete [] component;
            component = new T[v.dimension];
        }
        for(int i = 0; i < v.dimension; i++)
            component[i] = v.component[i];
        dimension = v.dimension;
    }
    return *this;
}  //  assignment operator

template<class T>
const dynamicVector<T>& dynamicVector<T>::operator=(const T& a){
    for(int i = 0; i < dimension; i++)
        component[i] = a;
    return *this;
}  //  assignment operator with a scalar argument

template<class T>
const dynamicVector<T>& dynamicVector<T>::operator+=( const dynamicVector<T>&v){
    for(int i = 0; i < dimension; i++)
        component[i] += v[i];
    return *this;
}  //  adding a dynamicVector to the current dynamicVector

template<class T>
const dynamicVector<T>& dynamicVector<T>::operator-=( const dynamicVector<T>&v){
    for(int i = 0; i < dimension; i++)
        component[i] -= v[i];
    return *this;
}  //  subtracting a dynamicVector from the current dynamicVector

template<class T>
const dynamicVector<T>& dynamicVector<T>::operator*=(const T& a){
    for(int i = 0; i < dimension; i++)
        component[i] *= a;
    return *this;
}  //  multiplying the current dynamicVector by a scalar


template<class T>
const dynamicVector<T>& dynamicVector<T>::operator/=(const T& a){
    for(int i = 0; i < dimension; i++)
        component[i] /= a;
    return *this;
}  //  dividing the current dynamicVector by a scalar

template<class T>
const dynamicVector<T> operator+(const dynamicVector<T>&u, const dynamicVector<T>&v){
    return dynamicVector<T>(u) += v;
}  //  dynamicVector plus dynamicVector

template<class T>
const dynamicVector<T> operator-(const dynamicVector<T>&u, const dynamicVector<T>&v){
    return dynamicVector<T>(u) -= v;
}  //  dynamicVector minus dynamicVector

template<class T>
const dynamicVector<T> operator*(const dynamicVector<T>&u, const T& a){
    return dynamicVector<T>(u) *= a;
}  //  dynamicVector times scalar

template<class T>
const dynamicVector<T> operator*(const T& a, const dynamicVector<T>&u){
    return dynamicVector<T>(u) *= a;
}  //  T times dynamicVector

template<class T>
const dynamicVector<T> operator*(double a, const dynamicVector<T>&u){
    return dynamicVector<T>(u) *= a;
}  //  scalar times dynamicVector

template<class T>
const dynamicVector<T> operator/(const dynamicVector<T>&u, const T& a){
    return dynamicVector<T>(u) /= a;
}  //  dynamicVector divided by scalar

template<class T>
const dynamicVector<T> operator-(const dynamicVector<T>&u){
    return dynamicVector<T>(u) *= -1.;
}  //  negative of a dynamicVector

template<class T>
T operator*(const dynamicVector<T>&u, const dynamicVector<T>&v){
    T sum = 0;
    for(int i = 0; i < u.dim(); i++)
        sum += u[i] * v[i];
    return sum;
}  //  dynamicVector times dynamicVector (inner product)

template<class T>
const dynamicVector<T>
operator|(const dynamicVector<T>&u, const dynamicVector<T>&v){
    dynamicVector<T> uv(u.dim(),0.);
    for(int i=0; i<u.dim(); i++)
        uv(i) = u[i] * v[i];
    return uv;
}  //  multiply component by component

template<class T>
ostream& operator <<(ostream& out, const dynamicVector<T>&v){
    //    for(int i = 0;i < v.dim(); i++)
    //        printf("v[%d]=%f;  ",i,(double)v[i]);
    //    printf("\n");
    out << "(";
    for (int i = 0; i < v.dim(); i++){
        if(i != v.dim() - 1)
            out << v[i] << ", ";
        else
            out << v[i] << ")";
    }
    return out;
}  //  print a dynamicVector

template<class T>
const T max(const dynamicVector<T>&u){
    T temp = u[0];
    for(int i = 0; i<u.dim(); i++)
        temp = temp >= u[i]? temp: u[i];
    return temp;
}  //  compute the max value of vector u

template<class T>
const T min(const dynamicVector<T>&u){
    T temp = u[0];
    for(int i = 0; i<u.dim(); i++)
        temp = temp <= u[i]? temp: u[i];
    return temp;
}  //  compute the min value of vector u

template<class T>
const T absMax(const dynamicVector<T>&u){
    T temp = fabs(u[0]);
    for(int i = 0; i<u.dim(); i++)
        temp = temp >= fabs(u[i]) ? temp: fabs(u[i]);
    return temp;
}  //  compute the max abs value of vector u

/*************************
 
 DFT
 
 ************************/
//template<class T, class S>
//const dynamicVector<S>
//FT(int inverse, const dynamicVector<T>&f, const S&type){
//    dynamicVector<S> result(f.dim(),0.);
//    double angle = 8. * atan(1.) / f.dim();
//    const complex w(cos(angle),(2*inverse-1) * sin(angle));
//    complex wi=1.;
//    polynomial<T> ff(f.dim());
//    for(int i=0; i<f.dim(); i++)
//        ff(i) = f[i];
//    for(int i=0; i<f.dim(); i++){
//        result(i) = ff(wi);
//        wi *= w;
//    }
//    return inverse ?
//    result / result.dim()
//    :
//    result;
//}  //  naive Fourier transform

/***************************
 
 FFT
 
 **************************/
template<class T, class S>
const dynamicVector<S>
FFT(int inverse, const dynamicVector<T>&f, const S&type, int first=1){
    if(f.dim()==1)return dynamicVector<S>(1,S(1.) * f[0]);
    dynamicVector<T> even(f.dim()/2,0.);
    dynamicVector<T> odd(f.dim()/2,0.);
    for(int i=0; i<f.dim(); i+=2){
        even(i/2) = f[i];
        odd(i/2) = f[i+1];
    }
    double angle = 8. * atan(1.) / f.dim();
    const complex w(cos(angle),(2*inverse-1) * sin(angle));
    complex wi = 1.;
    const dynamicVector<S> Feven = FFT(inverse,even,type, 0);
    dynamicVector<S> Fodd = FFT(inverse,odd,type, 0);
    for(int i=0; i<f.dim()/2; i++){
        Fodd(i) *= wi;
        wi *= w;
        //      Mults += 2;
    }
    return inverse&&first ?
    (1./f.dim()) * dynamicVector<S>(Feven + Fodd,Feven - Fodd)
    :
    dynamicVector<S>(Feven + Fodd,Feven - Fodd);
    
}  //  FFT

/***********************
 
 Convolution
 
 **********************/
template<class T>
const dynamicVector<complex>
conv(const dynamicVector<T>&u, const dynamicVector<T>&v){
    const complex I(0.,1.);
    dynamicVector<complex> uv(u.dim()+v.dim()-1,0.);
    int power2 = 1;
    for(; power2 < uv.dim(); power2 *= 2)
        ;
    dynamicVector<T> U(power2,0.);
    dynamicVector<T> V(power2,0.);
    for(int i=0; i<u.dim(); i++)
        U(i) = u[i];
    for(int i=0; i<v.dim(); i++)
        V(i) = v[i];
    dynamicVector<complex> UV = FFT(1,FFT(0,U,I) | FFT(0,V,I),I);
    for(int i = 0; i < uv.dim(); i++)
        uv(i) = UV[i];
    return uv;
}  //  convolution


#endif
