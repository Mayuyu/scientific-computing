//
//  main.cpp
//  scientific-computing
//
//  Created by 马 征 on 13-10-17.
//  Copyright (c) 2013年 马 征. All rights reserved.
//

#include <iostream>
#include "gauss_wgts.h"
#include <cmath>
using namespace std;

template <class T>
T harmonic_F(T q, T a, T r, T r_s, T eta, T epsilon_o, T epsilon_i) {
    double EPS=10e-4;
    T x,R,t,sum=0.0,p3,p2=0.0,p1=1.0,e=1.0;
    int i=0;
    x=cos(eta);
    t=R=(a*a)/(r*r_s);
    while (abs(e)>EPS) {
        p3=p2;
        p2=p1;
        p1=((2.0*i+1.0)*x*p2-i*p3)/(i+1);
        t=R*t;
        e=(t*(i+1)*(epsilon_o-epsilon_i)*p1)/((i+1)*epsilon_i+(i+2)*epsilon_o);
        sum+=e;
        i=i+1;
    }
    return (sum*q)/(a*epsilon_o);
}

template <class T>
dynamicVector<T> f(const dynamicVector<T> u, T r, T eta, T sigma) {
    dynamicVector<T> tmp(u.dim(),0.);
    for (int i =0; i<u.dim(); i++) {
        tmp(i)=1.0/sqrt(r*r+pow(u[i], 2.0/sigma)-2.0*r*pow(u[i], 1.0/sigma)*cos(eta));
    }
    return tmp;
}

template <class T>
T image_F(T q, T a, T r, T r_s, T eta, T epsilon_o, T epsilon_i) {
    T gamma, sigma, q_k, r_k, integral;
    gamma=(epsilon_i-epsilon_o)/(epsilon_i+epsilon_o);
    sigma=0.5*(1.0-gamma);
    r_k=a*a/r_s;
    q_k=-gamma*a/r_s;
    dynamicVector<T> w(1000,0.0), x(1000,0.0);
    gauleg(0.0, pow(r_k,sigma), x, w);
    integral=w*f(x,r, eta, sigma);
    return q_k/(epsilon_o*sqrt(r*r+r_k*r_k-2.0*r*r_k*cos(eta)))+integral*gamma*pow(r_k, 1.0-sigma)/(epsilon_o*a);
}

int main(int argc, const char * argv[])
{
    cout<<harmonic_F(1.0, 1.0, 2.0, 2.0, 0.0, 1.0, 80.0)<<endl;
    cout<<image_F(1.0, 1.0, 2.0, 2.0, 0.0, 1.0, 80.0)<<endl;
    return 0;
}

