//
//  main.cpp
//  scientific-computing
//
//  Created by 马 征 on 13-10-17.
//  Copyright (c) 2013年 马 征. All rights reserved.
//

#include <iostream>
#include <cmath>
#include "complex.h"
#include "dynamicMatrix.h"
#include "fftw3.h"
#include <fstream>

using namespace std;



/*************************************
 
    Project 4
 
 *******************************/


template <class T>
void split(const T& u, const T& k, const T& rho, dynamicMatrix<T>& ap, dynamicMatrix<T>& am) {
    dynamicMatrix<T> ip(2,2,0.0), p(2,2,0.0), eigp(2,2,0.0), eigm(2,2,0.0);
    T t1=sqrt(k/rho), t2=sqrt(k*rho);
    ip(0,0)=0.5; ip(0,1)=0.5; ip(1,0)=0.5/t2; ip(1,1)=-0.5/t2;
    p(0,0)=1.0; p(0,1)=t2; p(1,0)=1.0; p(1,1)=-t2;
    eigp(0,0)=(u+t1)>0 ? (u+t1):0.0; eigp(1,1)=(u-t1)>0 ? (u-t1):0.0;
    eigm(0,0)=(u+t1)<0 ? (u+t1):0.0; eigm(1,1)=(u-t1)<0 ? (u-t1):0.0;
    ap=ip*eigp*p;
    am=ip*eigm*p;
}



int main(int argc, const char * argv[]) {
    const int N=100, T=200;
    double u0=1.0, k=1.0, rho=1.0, t_delta=1.0/T, x_delta=1./N, r;
    r=t_delta/x_delta;
    dynamicMatrix<double> ap(2,2,0.0), am(2,2,0.0);
    dynamicVector<double> u(N+1,0.0), p(N+1,0.0), utmp(N+1,0.0), ptmp(N+1,0.0);
    for (int i=0; i<N+1; i++) {
        if (i*x_delta>0.4 && i*x_delta<0.6) {
            p(i)=1.0;
        }
    }
    split(u0, k, rho, ap, am);
    for (int t=0; t<T+1; t++) {
        for (int i=0; i<N+1; i++) {
            ptmp(i)=p(i)-r*(ap(0,0)*(p[i]-p[(i+N)%(N+1)])+ap(0,1)*(u[i]-u[(i+N)%(N+1)])+am(0,0)*(p[(i+1)%(N+1)]-p[i])+am(0,1)*(u[(i+1)%(N+1)]-u[i]));
            utmp(i)=u(i)-r*(ap(1,0)*(p[i]-p[(i+N)%(N+1)])+ap(1,1)*(u[i]-u[(i+N)%(N+1)])+am(1,0)*(p[(i+1)%(N+1)]-p[i])+am(1,1)*(u[(i+1)%(N+1)]-u[i]));
        }
        p=ptmp;
        u=utmp;
    }
    ofstream fout1("p4.txt");
    for (int i=0; i<N+1; i++) {
        fout1 << i*x_delta << " " << p[i] << endl;
    }
    fout1.close();
    cout << p << endl;
    return 0;
}






















///********************************************
// 
//            Project 3
// 
// *******************************************/
//
//double sq(double x) {
//    return x*x;
//}
//
//template <class T>
//void Init(const dynamicVector<double> &w, const dynamicVector<T> &a, dynamicVector<T> &tau,dynamicVector<double> &c, const double e, const int m) {
//    int N,uk,qq;
//    double b,P;
//    N=a.dim()-1;
////    b=log(1.0/e);
////    qq=2*b*M_PI;
//    qq=14;
//    b=1.5629;
////    qq=5;
////    b=0.5993;
//    for (int j=-N/2; j<N/2+1; j++) {
//        c(j+N/2)=exp(b*sq(2.0*M_PI*j/(m*N)));
//    }
//    for (int k=0; k<N+1; k++) {
//        uk=nearbyint(m*w[k]);
//        for (int j=-qq; j<qq+1; j++) {
//            P=0.5*sqrt(1.0/(b*M_PI))*exp(-sq(m*w[k]-(uk+j))/(4.0*b));
//            if (uk+j+m*N/2<0) {
//                tau((uk+j+m*N/2)%(m*N)+m*N)=tau[(uk+j+m*N/2)%(m*N)+m*N]+P*a[k];
//            }
//            else {
//                tau((uk+j+m*N/2)%(m*N))=tau[(uk+j+m*N/2)%(m*N)]+P*a[k];
//            }
//        }
//    }
//}
//
//template <class T>
//dynamicVector<T> convert(const dynamicVector<T> &x) {
//    int n=x.dim();
//    dynamicVector<T> tmp(n, 0);
//    for (int i=0; i<n; i++) {
//        if (i<n/2) {
//            tmp(i)=x[i+n/2];
//        } else {
//            tmp(i)=x[i-n/2];
//        }
//    }
//    return tmp;
//}
//
//template <class T>
//const double err_infty(const dynamicVector<T>& f, const dynamicVector<T>& a) {
//    int n=a.dim();
//    double sum=0.0;
//    for (int i=0; i<n; i++) {
//        sum+=sqrt(abs2(a[i]));
//    }
//    return absMax(f)/sum;
//}
//
//template <class T>
//const double err2(const dynamicVector<T>& f1, const dynamicVector<T>& f) {
//    int n=f.dim();
//    double sum1=0.0, sum2=0.0;
//    for (int i=0; i<n; i++) {
//        sum1+=abs2(f1[i]-f[i]);
//        sum2+=abs2(f[i]);
//    }
//    return sqrt(sum1/sum2);
//}
//
//int main(int argc, const char * argv[]) {
//    clock_t start, end, start1, end1,start2, end2, start3, end3;
//    const int N=4096, m=2;
//    dynamicVector<complex> a(N+1,0.),tau(m*N,0.),f(N+1,0.);
//    dynamicVector<double> w(N+1,0.),c(N+1,0.);
//    complex I(0.,1.);
//    srand((unsigned)time(NULL));
//    for (int i=0; i<N+1; i++) {
//        w(i)=(double)(rand()/(double)(RAND_MAX/N))-0.5*N;
//        complex random((double)(rand()/(double)RAND_MAX),(double)(rand()/(double)RAND_MAX));
//        a(i)=random;
//    }
//    
//    /**************************
//     
//        Direct DFT
//     
//     *************************/
//    
//    start=clock();
//    for (int i=0; i<N+1; i++) {
//        complex sum=0;
//        for (int j=0; j<N+1; j++) {
//            double tmp=w[j]*2.0*M_PI*(i-N/2)/N;
//            complex tmpc(cos(tmp),sin(tmp));
//            sum+=a[j]*tmpc;
//        }
//        f(i)=sum;
//    }
//    end=clock();
//    cout << "Direct method takes " << (double)(end-start)/CLOCKS_PER_SEC << " second." << endl;
//    
//    /*************************
//     
//        Initial
//     
//     ************************/
//    
//    start1=clock();
//    Init(w, a, tau, c,10e-6, m);
//    end1=clock();
//    cout << "Initial takes " << (double)(end1-start1)/CLOCKS_PER_SEC << " second." << endl;
//
//    /**********************
//     
//        FFT vesion 1
//     
//     *********************/
//    
//    dynamicVector<complex> f1(N+1,0.),tmp(m*N,0.), tmp1(m*N,0.), tmp2(m*N,0.);
//    tmp1=convert(tau);
//    start2=clock();
//    tmp2=FFT(1, tmp1, I);
//    end2=clock();
//    cout << "FFT vesion 1 takes " << (double)(end2-start2)/CLOCKS_PER_SEC << " second." << endl;
//    tmp=convert(m*N*tmp2);
//    for (int i=0; i<N+1; i++) {
//        f1(i)=c[i]*tmp[i+(m-1)*N/2];
//    }
//    cout << "The err_infty for FFT version 1: " << err_infty(f1-f, a) << endl;
//    cout << "The error2 for FFT version 1: " << err2(f1, f) << endl;
//    
//    /**********************
//     
//     FFT vesion FFTW
//     
//     *********************/
//    
//    dynamicVector<complex> tw(m*N,0.),f1w(m*N,0.),f2w(m*N,0.),fw(N+1,0.);
//    tw=convert(tau);
//    fftw_complex *in, *out;
//    fftw_plan p;
//    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(m*N));
//    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(m*N));
//    for (int i=0; i<m*N; i++) {
//        in[i][0]=tw[i].re();
//        in[i][1]=tw[i].im();
//    }
//    p = fftw_plan_dft_1d(m*N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
//    start3=clock();
//    fftw_execute(p); /* repeat as needed */
//    end3=clock();
//    cout << "FFTw takes " << (double)(end3-start3)/CLOCKS_PER_SEC << " second." << endl;
//    for (int i=0; i<m*N; i++) {
//        complex tmp(out[i][0],out[i][1]);
//        f1w(i)=tmp;
//    }
//    fftw_destroy_plan(p);
//    fftw_free(in); fftw_free(out);
//    f2w=convert(f1w);
//    for (int i=0; i<N+1; i++) {
//        fw(i)=c[i]*f2w[i+(m-1)*N/2];
//    }
//    cout << "The err_infty for FFTw: " << err_infty(fw-f, a) << endl;
//    cout << "The error2 for FFTw: " << err2(fw, f) << endl;
//    cout << err_infty(f1-f, a) << " " << err2(f1, f) << " " << err_infty(fw-f, a) << " " << err2(fw, f) << " " << (double)(end1-start1)/CLOCKS_PER_SEC << " " << (double)(end2-start2)/CLOCKS_PER_SEC << " " << (double)(end3-start3)/CLOCKS_PER_SEC << " " << (double)(end-start)/CLOCKS_PER_SEC << endl;
//    return 0;
//}


/**************************************************
 
            Project 2
 
 *************************************************/

//double pi =3.141592654;
//
//
//double sigma(double x, double y){
// //  return 1.0;
//    if (sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5))<0.25) {
//        return 0.1;
//    }
//    else
//        return 1.0;
//}
//
//double fxy(double x, double y){
////    return 2.0*pi*pi*sin(x*pi)*sin(y*pi);
//    return exp(-10.0*(x-0.1)*(x-0.1)-10.0*(y-0.1)*(y-0.1));
//}
//
//int main(int argc, const char * argv[])
//{
//    const int N=100;
//    
//    
//    dynamicVector<double> f((N-1)*(N-1), 0.0),x((N-1)*(N-1), 0.0), y((N-1)*(N-1), 0.0);
//    double delta=1.0/N;
//    sparseMat<double> u((N-1)*(N-1),(N-1)*(N-1), (5*N-9)*(N-1));
//    int t=0;
//    for (int k=0; k<(N-1)*(N-1); k++) {
//        int i=k/(N-1)+1;
//        int j=k%(N-1)+1;
//        f(k)=-fxy(i*delta,j*delta)/(N*N);
//        y(k)=sin(pi*i*delta)*sin(pi*j*delta);
//        double s1,s2,s3,s4;
//        s1=sigma((i-0.5)*delta,j*delta);
//        s2=sigma(i*delta,(j-0.5)*delta);
//        s3=sigma(i*delta,(j+0.5)*delta);
//        s4=sigma((i+0.5)*delta,j*delta);
//        bool b=false;
//        if (i!=1) {
//            u.val(t)=s1;
//            u.row_ind(t)=k-N+1;
//            u.col_ptr(k)=t;
//            b=true;
//            t++;
//        }
//        if (j!=1) {
//            u.val(t)=s2;
//            u.row_ind(t)=k-1;
//            if (b==false) {
//                u.col_ptr(k)=t;
//                b=true;
//            }
//            t++;
//        }
//        u.val(t)=-(s1+s2+s3+s4);
//        u.row_ind(t)=k;
//        if (b==false) {
//            u.col_ptr(k)=t;
//            b=true;
//        }
//        t++;
//        if (j!=N-1) {
//            u.val(t)=s3;
//            u.row_ind(t)=k+1;
//            if (b==false) {
//                u.col_ptr(k)=t;
//                b=true;
//            }
//            t++;
//        }
//        if (i!=N-1) {
//            u.val(t)=s4;
//            u.row_ind(t)=k+N-1;
//            if (b==false) {
//                u.col_ptr(k)=t;
//                b=true;
//            }
//            t++;
//        }
//    }
//    u.col_ptr((N-1)*(N-1))=(5*N-9)*(N-1);
//
//    int iter=0;
//    double err=0.0;
//    sparseLinbcg<double> U(u);
//    time_t start,finish;
//    double duration;
//    start=clock();
//    U.solve(f, x, 1, 10e-5, 10000, iter, err);
//    finish=clock();
//    duration=(double)(finish - start) / CLOCKS_PER_SEC;
//    cout << iter << endl;
//    cout << duration << endl;
//    cout << err << endl;
//    ofstream fout01("bcg2_0.2_100.txt");
//    for (int i=0; i<N+1; i++) {
//        if (i==0 || i==N ) {
//            for (int j=0; j<N+1; j++) {
//                fout01 << i*delta << " " << j*delta << " " << 0 << endl;
//            }
//        } else {
//            for (int j=0; j<N+1; j++) {
//                if (j==0 || j==N) {
//                    fout01 << i*delta << " " << j*delta << " " << 0 << endl;
//                }
//                else
//                fout01 << i*delta << " " << j*delta << " " << x[(i-1)*(N-1)+j-1] << endl;
//                }
//            }
//    }
//    fout01.close();
    
        
    

/******************************************
    
                SOR
    
******************************************/
    
//    dynamicMatrix<double> a(N+1, N+1, 0.0),b(N+1, N+1, 0.0),c(N+1, N+1, 0.0),d(N+1, N+1, 0.0),e(N+1, N+1, 0.0),f1(N+1, N+1, 0.0),x1(N+1, N+1, 0.0);
//    for (int i=1; i<N; i++) {
//        for (int j=1; j<N; j++) {
//            double s1,s2,s3,s4;
//            s1=sigma((i+0.5)*delta,j*delta);
//            s2=sigma((i-0.5)*delta,j*delta);
//            s3=sigma(i*delta,(j+0.5)*delta);
//            s4=sigma(i*delta,(j-0.5)*delta);
//
//            a(i,j)=s1;
//            b(i,j)=s2;
//            c(i,j)=s3;
//            d(i,j)=s4;
//            e(i,j)=-(s1+s2+s3+s4);
//            f1(i,j)=-fxy(i*delta,j*delta)/(N*N);
//        }
//    }
//    
//    sor(a, b, c, d, e, f1, x1, cos(pi/double(N)));
//    
//    cout << x1 << endl;
//    
//    ofstream fout02("sor2_0.1_50.txt");
//    for (int i=0; i<N+1; i++) {
//            for (int j=0; j<N+1; j++) {
//                fout02 << i*delta << " " << j*delta << " " << x1(i,j,"read") << endl;
//            }
//        }
//    
//    fout02.close();
//    
//    
//    return 0; 
//
//}




/**************************************************
 
                Project 1
 
 *************************************************/

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

