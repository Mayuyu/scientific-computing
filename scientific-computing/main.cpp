//
//  main.cpp
//  scientific-computing
//
//  Created by 马 征 on 13-10-17.
//  Copyright (c) 2013年 马 征. All rights reserved.
//

#include <iostream>
#include <fstream>
#include "gauss_wgts.h"
#include "linbcg.h"
#include "asolve.h"
#include "sparse.h"
#include <cmath>
#include "sor.h"
using namespace std;

double pi =3.141592654;


double sigma(double x, double y){
    return 1.0;
}

double fxy(double x, double y){
    return 2.0*pi*pi*sin(x*pi)*sin(y*pi);
}

int main(int argc, const char * argv[])
{
    const int N=4;
    dynamicVector<double> f((N-1)*(N-1), 0.0),x((N-1)*(N-1), 1.0), y((N-1)*(N-1), 0.0);
    double delta=1.0/N;
    sparseMat<double> u((N-1)*(N-1),(N-1)*(N-1), (5*N-9)*(N-1));
    int t=0;
    for (int k=0; k<(N-1)*(N-1); k++) {
        int i=k/(N-1)+1;
        int j=k%(N-1)+1;
        f(k)=-fxy(i*delta,j*delta)/(N*N);
        y(k)=sin(pi*i*delta)*sin(pi*j*delta);
        double s1,s2,s3,s4;
        s1=sigma((i-0.5)*delta,j*delta);
        s2=sigma(i*delta,(j-0.5)*delta);
        s3=sigma(i*delta,(j+0.5)*delta);
        s4=sigma((i+0.5)*delta,j*delta);
        bool b=false;
        if (i!=1) {
            u.val(t)=s1;
            u.row_ind(t)=k-N+1;
            u.col_ptr(k)=t;
            b=true;
            t++;
        }
        if (j!=1) {
            u.val(t)=s2;
            u.row_ind(t)=k-1;
            if (b==false) {
                u.col_ptr(k)=t;
                b=true;
            }
            t++;
        }
        u.val(t)=-(s1+s2+s3+s4);
        u.row_ind(t)=k;
        if (b==false) {
            u.col_ptr(k)=t;
            b=true;
        }
        t++;
        if (j!=N-1) {
            u.val(t)=s3;
            u.row_ind(t)=k+1;
            if (b==false) {
                u.col_ptr(k)=t;
                b=true;
            }
            t++;
        }
        if (i!=N-1) {
            u.val(t)=s4;
            u.row_ind(t)=k+N-1;
            if (b==false) {
                u.col_ptr(k)=t;
                b=true;
            }
            t++;
        }
    }
    u.col_ptr((N-1)*(N-1))=(5*N-9)*(N-1);

    int iter=0;
    double err=0.0;
    sparseLinbcg<double> U(u);
    U.solve(f, x, 1, 1.0e-6, 100, iter, err);
cout << f <<endl;
    cout << x <<endl;
 //   cout << u.val << endl;
 //   cout << t << endl;
//    cout << delta << endl;
    cout << y << endl;
 //   cout << u.row_ind << endl;
 //   cout << u.col_ptr << endl;
    cout << x-y << endl;
    cout << err << endl;
    cout << iter << endl;
    cout << 2.0*pi*pi/4.0 << endl;
    
/******************************************
    
                SOR
    
******************************************/
    
    dynamicMatrix<double> a(N+1, N+1, 0.0),b(N+1, N+1, 0.0),c(N+1, N+1, 0.0),d(N+1, N+1, 0.0),e(N+1, N+1, 0.0),f1(N+1, N+1, 0.0),x1(N+1, N+1, 0.0);
    for (int i=1; i<N; i++) {
        for (int j=1; j<N; j++) {
            double s1,s2,s3,s4;
            s1=sigma((i+0.5)*delta,j*delta);
            s2=sigma((i-0.5)*delta,j*delta);
            s3=sigma(i*delta,(j+0.5)*delta);
            s4=sigma(i*delta,(j-0.5)*delta);

            a(i,j)=s1;
            b(i,j)=s2;
            c(i,j)=s3;
            d(i,j)=s4;
            e(i,j)=-(s1+s2+s3+s4);
            f1(i,j)=-fxy(i*delta,j*delta)/(N*N);
        }
    }
    
    sor(a, b, c, d, e, f1, x1, cos(pi/double(N)));
    
    cout << x1 << endl;
    
    
    
    return 0;

}





















//template <class T>
//T harmonic_F(T q, T a, T r, T r_s, T eta, T epsilon_o, T epsilon_i, int &i) {
//    double EPS=10e-8;
//    T x,R,t,sum=0.0,p3,p2=0.0,p1=1.0,e=1.0;
//    i=0;
//    x=cos(eta);
//    t=R=(a*a)/(r*r_s);
////    while (abs(e)>EPS) {
////        p3=p2;
////        p2=p1;
////        p1=((2.0*i+1.0)*x*p2-i*p3)/(i+1);
////        t=R*t;
////        e=(t*(i+1)*(epsilon_o-epsilon_i)*p1)/((i+1)*epsilon_i+(i+2)*epsilon_o);
////        sum+=e;
////        i=i+1;
////    }
//    while (i<1000) {
//        p3=p2;
//        p2=p1;
//        p1=((2.0*i+1.0)*x*p2-i*p3)/(i+1);
//        t=R*t;
//        e=(t*(i+1)*(epsilon_o-epsilon_i)*p1)/((i+1)*epsilon_i+(i+2)*epsilon_o);
//        sum+=e;
//        i=i+1;
//    }
//
//    return (sum*q)/(a*epsilon_o);
//}
//
//template <class T>
//dynamicVector<T> f(const dynamicVector<T> u, T r, T r_k, T eta) {
//    dynamicVector<T> tmp(u.dim(),0.);
//    for (int i =0; i<u.dim(); i++) {
//        tmp(i)=1.0/sqrt(r*r+0.25*r_k*r_k*(1.0+u[i])*(1.0+u[i])-r_k*r*(1.0+u[i])*cos(eta));
//    }
//    return tmp;
//}
//
//template <class T>
//T image_F(T q, T a, T r, T r_s, T eta, T epsilon_o, T epsilon_i) {
//    double EPS=10e-6;
//    T gamma, sigma, q_k, r_k, integral=0,tmp=0;
//    int i=10;
//    gamma=(epsilon_i-epsilon_o)/(epsilon_i+epsilon_o);
//    sigma=0.5*(1.0-gamma);
//    r_k=a*a/r_s;
//    q_k=-gamma*a/r_s;
////    do {
////        tmp=integral;
////        dynamicVector<T> w(i,0.0), x(i,0.0);
////        gaujac(x,w,0.0,sigma-1.0);
////        integral=w*f(x,r,r_k,eta);
////        i++;
////    } while (abs(integral-tmp)>EPS);
//    dynamicVector<T> w(i,0.0), x(i,0.0);
//    gaujac(x,w,0.0,sigma-1.0);
//    integral=w*f(x,r,r_k,eta);
//    return q_k/(epsilon_o*sqrt(r*r+r_k*r_k-2.0*r*r_k*cos(eta)))+integral*gamma*sigma*r_k/(epsilon_o*a*pow(2.0, sigma));
//}
//
////double har_F(double r_a, double r_b, double eta) {
////    return harmonic_F(1.0, 1.0, r_a, r_b, eta, 1.0, 80.0);
////}
////
////double ima_F(double r_a, double r_b, double eta) {
////    return image_F(1.0, 1.0, r_a, r_b, eta, 1.0, 80.0);
////}
////double har_G(double r_a, double r_b, double eta) {
////    return harmonic_F(1.0, 1.0, r_a, r_b, eta, 1.0, 80.0)+1.0/(sqrt(r_a*r_a+r_b*r_b-2.0*r_a*r_b*cos(eta)));
////}
////
////double ima_G(double r_a, double r_b, double eta) {
////    return image_F(1.0, 1.0, r_a, r_b, eta, 1.0, 80.0)+1.0/(sqrt(r_a*r_a+r_b*r_b-2.0*r_a*r_b*cos(eta)));
////}
//int main(int argc, const char * argv[])
//{
////    ofstream fout01("Us_har.txt");
////    for (int i=1; i<100; i++) {
////        fout01 << 1.0+i*0.005 << " " << 0.5*har_F(1.0+i*0.005, 1.0+i*0.005, 0.0) << endl;
////    }
////    fout01.close();
////    ofstream fout02("Us_ima.txt");
////    for (int i=1; i<100; i++) {
////        fout02 << 1.0+i*0.005 << " " << 0.5*ima_F(1.0+i*0.005, 1.0+i*0.005, 0.0) << endl;
////    }
////    fout02.close();
////    ofstream fout03("Up_har0.txt");
////    for (int i=1; i<100; i++) {
////        fout03 << 1.0+i*0.005 << " " << har_G(1.0+i*0.005, 1.0, 0.0) << endl;
////    }
////    ofstream fout04("Up_harpi.txt");
////    for (int i=1; i<100; i++) {
////        fout04 << 1.0+i*0.005 << " " << har_G(1.0+i*0.005, 1.0, M_PI) << endl;
////    }
////    fout04.close();
////    ofstream fout05("Up_ima0.txt");
////    for (int i=1; i<100; i++) {
////        fout05 << 1.0+i*0.005 << " " << ima_G(1.0+i*0.005, 1.0, 0.0) << endl;
////    }
////    fout05.close();
////    ofstream fout06("Up_imapi.txt");
////    for (int i=1; i<100; i++) {
////        fout06 << 1.0+i*0.005 << " " << ima_G(1.0+i*0.005, 1.0, M_PI) << endl;
////    }
////    fout06.close();
////    ofstream fout07("Up_har2.txt");
////    for (int i=0; i<361; i++) {
////        fout07 << i*M_PI/180.0 << " " << har_G(2.0, 1.0, i*M_PI/180.0) << endl;
////    }
////    fout07.close();
////    ofstream fout08("Up_ima2.txt");
////    for (int i=0; i<361; i++) {
////        fout08 << i*M_PI/180.0 << " " << ima_G(2.0, 1.0, i*M_PI/180.0) << endl;
////    }
////    fout08.close();
//    int n;
//    ofstream fout09("truc_num.txt");
//    for (int i=1; i<100; i++) {
//        harmonic_F(1.0, 1.0, 1.0+i*0.005, 1.0+i*0.005, 0.0, 1.0, 80.0,n);
//        fout09 << 1.0+i*0.005 << " " << n << endl;
//    }
//    fout09.close();
//    return 0;
//}

