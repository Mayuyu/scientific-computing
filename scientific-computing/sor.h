//
//  sor.h
//  scientific-computing
//
//  Created by 马 征 on 13-11-6.
//  Copyright (c) 2013年 马 征. All rights reserved.
//

#ifndef scientific_computing_sor_h
#define scientific_computing_sor_h

#include "dynamicMatrix.h"

template <class T>
void sor(const dynamicMatrix<T> &a, const dynamicMatrix<T> &b, const dynamicMatrix<T> &c, const dynamicMatrix<T> &d, const dynamicMatrix<T> &e, const dynamicMatrix<T> &f, dynamicMatrix<T> &u, const double rjac)
//Successive overrelaxation solution of equation (20.5.25) with Chebyshev acceleration. a, b, c, d, e, and f are input as the coefficients of the equation, each dimensioned to the grid size [0..jmax-1][0..jmax-1]. u is input as the initial guess to the solution, usually zero, and returns with the final value. rjac is input as the spectral radius of the Jacobi iteration, or an estimate of it.
{
    const int MAXITS=1000;
    const double EPS=1.0e-13;
    double anormf=0.0,omega=1.0;
    int jmax=a.height();
    for (int j=1;j<jmax-1;j++) //Compute initial norm of residual and terminate iterations when norm has been reduced by a factor EPS.
        for (int l=1;l<jmax-1;l++)
            anormf += abs(f(j,l, "read"));
            for (int n=0;n<MAXITS;n++) {
                double anorm=0.0;
                int jsw=1;
                for (int ipass=0;ipass<2;ipass++) {
                    int lsw=jsw;
                    for (int j=1;j<jmax-1;j++) {
                       // Assumes initial u is zero.Odd-even ordering.
                        for (int l=lsw;l<jmax-1;l+=2) {
                            double resid=a(j,l,"read")*u(j+1,l,"read")+b(j,l,"read")*u(j-1,l,"read")+c(j,l,"read")*u(j,l+1,"read")+d(j,l,"read")*u(j,l-1,"read")+e(j,l,"read")*u(j,l,"read")-f(j,l,"read");
                            anorm += abs(resid);
                            u(j,l) -= omega*resid/e(j,l,"read");
                        }
                        lsw=3-lsw;
                    }
                    jsw=3-jsw;
                    omega=(n == 0 && ipass == 0 ? 1.0/(1.0-0.5*rjac*rjac) :
                           1.0/(1.0-0.25*rjac*rjac*omega));
                    }
                    if (anorm < EPS*anormf) return;
                }
                throw("MAXITS exceeded");
}

#endif
