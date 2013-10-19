//
//  gauss_wgts.h
//  scientific-computing
//
//  Created by 马 征 on 13-10-17.
//  Copyright (c) 2013年 马 征. All rights reserved.
//

#ifndef scientific_computing_gauss_wgts_h
#define scientific_computing_gauss_wgts_h

#include "dynamicVector.h"
#include "gamma.h"

template <class T>
void gauleg(const T x1, const T x2, dynamicVector<T> &x, dynamicVector<T> &w)
//Given the lower and upper limits of integration x1 and x2, this routine returns arrays x[0..n-1] and w[0..n-1] of length n, containing the abscissas and weights of the Gauss-Legendre n-point quadrature formula.
{
    const double EPS=1.0e-14;                //EPS is the relative precision
    T z1,z,xm,xl,pp,p3,p2,p1;
    int n=x.dim();
    int m=(n+1)/2;
    xm=0.5*(x2+x1);                          //The roots are symmetric in the interval, so we only have to find half of them.
    xl=0.5*(x2-x1);
    for (int i=0;i<m;i++) {                  //Loop over the desired roots.
        z=cos(3.141592654*(i+0.75)/(n+0.5)); //Starting with this approximation to the ith root, we enter the main loop of refinement by Newton’s method.
        do {
            p1=1.0;
            p2=0.0;
            for (int j=0;j<n;j++) {           //Loop up the recurrence relation to get the Legendre polynomial evaluated at z.
                p3=p2;
                p2=p1;
                p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
            }
            //p1 is now the desired Legendre polynomial. We next compute pp, its derivative, by a standard relation involving also p2, the polynomial of one lower order.
            pp=n*(z*p1-p2)/(z*z-1.0);
            z1=z;                             //Newton’s method.
            z=z1-p1/pp;
        } while (abs(z-z1) > EPS);
        x(i)=xm-xl*z;                         //Scale the root to the desired interval, and put in its symmetric counterpart.
        x(n-1-i)=xm+xl*z;
        w(i)=2.0*xl/((1.0-z*z)*pp*pp);        //Compute the weight
        w(n-1-i)=w[i];                        //and its symmetric conterpart
    }
}

template <class T>
void gaulag(dynamicVector<T> &x, dynamicVector<T> &w, const T alf)
//Given alf, the parameter  ̨ of the Laguerre polynomials, this routine returns arrays x[0..n-1] and w[0..n-1] containing the abscissas and weights of the n-point Gauss-Laguerre quadrature formula. The smallest abscissa is returned in x[0], the largest in x[n-1].
{
    const int MAXIT=10;                             //EPS is the relative precision.
    const double EPS=1.0e-14;
    int i,its,j;
    double ai,p1,p2,p3,pp,z,z1;
    int n=x.dim();
    for (i=0;i<n;i++) {
        if (i == 0) {                               //Loop over the desired roots.
            z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
        }
        else if(i == 1) {
            z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
        }
        else {
            ai=i-1;
            z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/(1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
        }
        for (its=0;its<MAXIT;its++) {               //Refinement by Newton’s method.
            p1=1.0;
            p2=0.0;
            for (j=0;j<n;j++) {                     //Look up the recurrence relation to get the Laguerre polynomial
                p3=p2;                              //evaluated at z
                p2=p1;
                p1=((2*j+1+alf-z)*p2-(j+alf)*p3)/(j+1);
            }
            pp=(n*p1-(n+alf)*p2)/z;
            z1=z;
            z=z1-p1/pp;                             //Newton’s formula.
            if (abs(z-z1) <= EPS) break;
        }
        if (its >= MAXIT) throw("too many iterations in gaulag");
        x(i)=z;                                     //Store the root and the weight.
        w(i) = -exp(gammln(alf+n)-gammln(double(n)))/(pp*n*p2);
    }
}

template <class T>
void gauher(dynamicVector<T> &x, dynamicVector<T> &w) //This routine returns arrays x[0..n-1] and w[0..n-1] containing the abscissas and weights of the n-point Gauss-Hermite quadrature formula. The largest abscissa is returned in x[0], the most negative in x[n-1].
{
    const double EPS=1.0e-14,PIM4=0.7511255444649425;
    const int MAXIT=10;                             //Maximum iterations.
    int i,its,j,m;
    double p1,p2,p3,pp,z,z1;
    int n=x.dim();
    m=(n+1)/2;                              //The roots are symmetric about the origin, so we have to find only half of them.
    for (i=0; i<m; i++) {
        if (i==0) {
             z=sqrt(double(2*n+1))-1.85575*pow(double(2*n+1),-0.16667);
        } else if (i==1) {
            z -= 1.14*pow(double(n),0.426)/z;
        } else if (i==2) {
            z=1.86*z-0.86*x[0];
        } else if (i==3) {
            z=1.91*z-0.91*x[1];
        } else {
            z=2.0*z-x[i-2];
        }
        for (its=0;its<MAXIT;its++) {
            p1=PIM4;
            p2=0.0;
            for (j=0;j<n;j++) {
                p3=p2;
                p2=p1;
                p1=z*sqrt(2.0/(j+1))*p2-sqrt(double(j)/(j+1))*p3;
            }
            pp=sqrt(double(2*n))*p2;
            z1=z;
            z=z1-p1/pp;                         //Newton’s formula.
            if (abs(z-z1) <= EPS) break;
        }
        if (its >= MAXIT) throw("too many iterations in gauher");
        x(i)=z;
        x(n-1-i)=-z;
        w(i)=2.0/(pp*pp);
        w(n-1-i)=w[i];
    }
}















#endif
