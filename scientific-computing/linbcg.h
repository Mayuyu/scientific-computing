//
//  linbcg.h
//  scientific-computing
//
//  Created by 马 征 on 13-11-5.
//  Copyright (c) 2013年 马 征. All rights reserved.
//

#ifndef scientific_computing_linbcg_h
#define scientific_computing_linbcg_h

#include "dynamicVector.h"
#include "sparse.h"

template<class T>
struct Linbcg {
    virtual void asolve(dynamicVector<T> &b, dynamicVector<T> &x, const int itrnsp) = 0;
    virtual void atimes(dynamicVector<T> &x, dynamicVector<T> &r, const int itrnsp) = 0;
    void solve(dynamicVector<T> &b, dynamicVector<T> &x, const int itol, const double tol,const int itmax, int &iter, double &err);
    double snrm(dynamicVector<T> &sx, const int itol);   //Utility used by solve.
};

template<class T>
void Linbcg<T>::solve(dynamicVector<T> &b, dynamicVector<T> &x, const int itol, const double tol,const int itmax, int &iter, double &err)
{
    double ak,akden,bk,bkden=1.0,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm;
    const double EPS=1.0e-14;
    int j,n=b.dim();
    dynamicVector<T> p(n),pp(n),r(n),rr(n),z(n),zz(n);
    iter=0;
    atimes(x,r,0);
    for (j=0;j<n;j++) {
        r(j)=b[j]-r[j];
        rr(j)=r[j];
    }
    
    //atimes(r,rr,0);      //Uncomment this line to get the “minimum residual” variant of the algorithm.
    if (itol == 1) {
        bnrm=snrm(b,itol);
        asolve(r,z,0);
    }
    else if (itol == 2) {
        asolve(b,z,0);
        bnrm=snrm(z,itol);
        asolve(r,z,0);
    }
    else if (itol == 3 || itol == 4) {
        asolve(b,z,0);
        bnrm=snrm(z,itol);
        asolve(r,z,0);
        znrm=snrm(z,itol);
    } else throw("illegal itol in linbcg");
    while (iter < itmax) {
        ++iter;
        asolve(rr,zz,1);
        for (bknum=0.0,j=0;j<n;j++) bknum += z[j]*rr[j];
        if (iter == 1) {
            for (j=0;j<n;j++) {
                p(j)=z[j];
                pp(j)=zz[j];
            }
        } else {
            bk=bknum/bkden;
            for (j=0;j<n;j++) {
                p(j)=bk*p[j]+z[j];
                pp(j)=bk*pp[j]+zz[j];
            }
        }
        bkden=bknum;
        atimes(p,z,0);
        for (akden=0.0,j=0;j<n;j++) akden += z[j]*pp[j];
        ak=bknum/akden;
        atimes(pp,zz,1);
        for (j=0;j<n;j++) {
            x(j) += ak*p[j];
            r(j) -= ak*z[j];
            rr(j) -= ak*zz[j];
        }
        asolve(r,z,0);
        if (itol == 1)
            err=snrm(r,itol)/bnrm;
        else if (itol == 2)
            err=snrm(z,itol)/bnrm;
        else if (itol == 3 || itol == 4) {
            zm1nrm=znrm;
            znrm=snrm(z,itol);
            if (abs(zm1nrm-znrm) > EPS*znrm) {
                dxnrm=abs(ak)*snrm(p,itol);
                err=znrm/abs(zm1nrm-znrm)*dxnrm;
            } else {
                err=znrm/bnrm;
                continue;
            }
            xnrm=snrm(x,itol);
            if (err <= 0.5*xnrm) err /= xnrm;
            else {
                err=znrm/bnrm;
                continue;
            }
        }
        if (err <= tol) break;
    }
}

template<class T>
double Linbcg<T>::snrm(dynamicVector<T> &sx, const int itol) {
    int i,isamax,n=sx.dim();
    double ans;
    if (itol <= 3) {
        ans = 0.0;
        for (i=0;i<n;i++) ans += sx[i]*sx[i];
        return sqrt(ans);
    } else {
        isamax=0;
        for (i=0;i<n;i++) {
            if (abs(sx[i]) > abs(sx[isamax])) isamax=i;
        }
        return abs(sx[isamax]);
    }
}







#endif
