//
//  3Dpoint.h
//  scientific-computing
//
//  Created by 马 征 on 14-1-7.
//  Copyright (c) 2014年 马 征. All rights reserved.
//

#ifndef scientific_computing__Dpoint_h
#define scientific_computing__Dpoint_h

#include <iostream>

using namespace std;

template<class T> class point{
    T x;
    T y;
    T z;
public:
    point(T x = 0, T y = 0, T z=0):x(x),y(y),z(z){}
    point(const point & u):x(u.x),y(u.y),z(u.z){}
    const T rx()const{return x;}
    const T ry()const{return y;}
    const T rz()const{return z;}
    const point& operator=(const point&);
    const point& operator+=(const point&);
    const point& operator-=(const point&);
    const point& operator*=(const point&);
    const point& operator/=(const T&);
    const point& operator/=(const point&);
    ~point(){}
};

template <class T>
const point<T>& point<T>::operator=(const point<T>& c){
    x = c.x;
    y = c.y;
    z = c.z;
    return *this;
}

template <class T>
const point<T>& point<T>::operator+=(const point<T>& c){
    x += c.x;
    y += c.y;
    z += c.z;
    return *this;
}

template <class T>
const point<T>& point<T>::operator-=(const point<T>& c){
    x -= c.x;
    y -= c.y;
    z -= c.z;
    return *this;
}
template <class T>
const point<T>& point<T>::operator/=(const T& d){
    x /= d;
    y /= d;
    z /= d;
    return *this;
}

template <class T>
T abs(const point<T>& c){
    return sqrt(c.rx()*c.rx()+c.ry()*c.ry()+c.rz()*c.rz());
}

template <class T>
const point<T> operator+(const point<T>& a, const point<T>& b){
    return point<T>(a)+=b;
}

template <class T>
const point<T> operator-(const point<T>& a, const point<T>& b){
    return point<T>(a)-=b;
}

template <class T>
const point<T> operator*(const T& a, const point<T>& b){
    point<T> tmp (a*b.rx(),a*b.ry(),a*b.rz());
    return tmp;
}

template <class T>
ostream& operator <<(ostream& out, const point<T>& c){
    out << "(" << c.rx() << "," << c.ry() << "," << c.rz() << ")";
    return out;
}

#endif
