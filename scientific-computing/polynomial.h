//
//  polynomial.h
//  scientific-computing
//
//  Created by 马 征 on 13-10-17.
//  Copyright (c) 2013年 马 征. All rights reserved.
//

#ifndef scientific_computing_polynomial_h
#define scientific_computing_polynomial_h

template<class T> class polynomial:public list<T>{
public:
    polynomial(int n=0){
        this->number = n;
        this->item = n ? new T*[n] : 0;
        for(int i=0; i<n; i++)
            this->item[i] = 0;
    }  //  constructor
    polynomial(int n, const T&a){
        this->number = n;
        this->item = n ? new T*[n] : 0;
        for(int i=0; i<n; i++)
            this->item[i] = new T(a);
    }  //  constructor with 'T' argument
    ~polynomial(){}
    int degree() const{return this->number-1;}
};

template<class T>
const T
calculatePolynomial(const polynomial<T>&p, const T&x){
    T powerOfX = 1;
    T sum=0;
    for(int i=0; i<p.size(); i++){
        sum += p[i] * powerOfX;
        powerOfX *= x;
    }
    return sum;
}  //  calculate a polynomial

template<class T>
const T
HornerPolynomial(const polynomial<T>&p, const T&x){
    T result = p[p.degree()];
    for(int i=p.degree(); i>0; i--){
        result *= x;
        result += p[i-1];
    }
    return result;
}  //  Horner algorithm to calculate a polynomial

template<class T>
const polynomial<T>& operator+=(polynomial<T>& p, const polynomial<T>&q){
    if(p.size() >= q.size())
        for(int i=0; i<q.size(); i++)
            p(i) += q[i];
    else{
        polynomial<T> keepQ = q;
        p = keepQ += p;
    }
    return p;
}  //  add polynomial

template<class T>
const polynomial<T> operator+(const polynomial<T>& p, const polynomial<T>&q){
    polynomial<T> keep = p;
    return keep += q;
}  //  add two polynomials

template<class T>
const polynomial<T>& operator*=(polynomial<T>& p, const T&a){
    for(int i=0; i<p.size(); i++)
        p(i) *= a;
    return p;
}  //  multiplication by scalar

template<class T>
const polynomial<T> operator*(const T&a, const polynomial<T>&p){
    polynomial<T> keep = p;
    return keep *= a;
}  //  scalar times polynomial

template<class T>
polynomial<T>
operator*(const polynomial<T>&p, const polynomial<T>&q){
    polynomial<T> result(p.degree()+q.degree()+1,0);
    for(int i=0; i<result.size(); i++)
        for(int j=max(0,i-q.degree());
            j<=min(i,p.degree()); j++){
            if(j == max(0,i-q.degree()))
                result(i) = p[j] * q[i-j];
            else
                result(i) += p[j] * q[i-j];
        }
    return result;
}  //  multiply 2 polynomials

template<class T>
polynomial<T>&
operator*=(polynomial<T>&p, const polynomial<T>&q){
    return p = p * q;
}  //  multiply by polynomial

template<class T>
const T
integral(const polynomial<T>&p){
    T sum = 0;
    for(int i=0; i<p.size(); i++)
        sum += (1./(i+1)) * p[i];
    return sum;
}  //  integral on the unit interval

template<class T>
const T
integral(const polynomial<polynomial<T> >&p){
    polynomial<T> sum(p.size()+1,0);
    polynomial<T> one(1,1);
    polynomial<T> x(2,0);
    x(1) = 1;
    polynomial<T> oneMinusX(2,1);
    oneMinusX(1) = -1;
    list<polynomial<T> > xPowers(p.size(),one);
    list<polynomial<T> > oneMinusXpowers(p.size()+1,one);
    for(int i=1; i<p.size(); i++)
        xPowers(i) = x * xPowers[i-1];
    for(int i=1; i<=p.size(); i++)
        oneMinusXpowers(i) = oneMinusX * oneMinusXpowers[i-1];
    for(int k=p.degree(); k>=0; k--)
        for(int j=0; j<=k; j++)
            sum += (p[k][j]/(j+1))
            * oneMinusXpowers[j+1] * xPowers[k-j];
    return integral(sum);
}  //  integral on the truangle

/***************************************
 
 Newtonian interpolation template
 
 ***************************************/
template<class T> class newtonian:public list<dynamicVector<T> >{
public:
    newtonian(const dynamicVector<T>&, const dynamicVector<T>&);
    const T& get(int i) const{return (*this->item[i])(i+1);}
    const T calculate(const T& x) const;
    int degree() const{return this->number-1;}
};

template<class T>
newtonian<T>::newtonian(const dynamicVector<T> &x, const dynamicVector<T> &y){
    this->number = x.dim();
    this->item = new dynamicVector<T>*[x.dim()];
    for (int i = 0; i < x.dim(); i++) {
        this->item[i] = new dynamicVector<T>(i + 2, x[i]);
    }
    for (int i = 1; i < x.dim()+1; i++) {
        for (int j = i-1; j < x.dim(); j++) {
            if (i==1) {
                (*this->item[j])(i) = y[j];
            }
            else{
                (*this->item[j])(i) = ((*this->item[j])(i-1)-(*this->item[j-1])(i-1))/((*this->item[j])(0)-(*this->item[j-i+1])(0));
            }
        }
    }
}

template<class T>
const T newtonian<T>::calculate(const T& x) const{
    T u = (*this->item[degree()])(degree()+1);
    for (int i = degree(); i > 0; i--) {
        u = (*this->item[i-1])(i) + (x - (*this->item[i-1])(0))*u;
    }
    return u;
}

template<class T>
const dynamicVector<T> calculate(const newtonian<T>& p, const dynamicVector<T>& u){
    dynamicVector<T> v(u);
    for (int i = 0; i < u.dim(); i++){
        v(i) = p.calculate(u[i]);
    }
    return v;
}

/**************************************************
 
 Cubic spline template
 
 ***************************************************/
template<class T> class spline:public list<dynamicVector<T> >{
public:
    spline(dynamicVector<T>&, dynamicVector<T>&);
    const T& get(int i) const{return (*this->item[i])(i+1);}
    const T calculate(const T& x) const;
    int degree() const{return this->number-1;}
};

template<class T>
spline<T>::spline(dynamicVector<T>& x, dynamicVector<T>& y){
    this->number = 3;
    this->item = new dynamicVector<T>*[3];
    this->item[0] = new dynamicVector<T>(x);
    this->item[1] = new dynamicVector<T>(y);
    this->item[2] = new dynamicVector<T>(x);
    dynamicVector<T> lamda(x.dim(),0.), miu(x.dim(),0.),beta(x.dim(),0.),alpha(x.dim(),0.);
    lamda(0) = 1.;
    miu(0) = 3.*(y(1)-y(0))/(x(1) -x(0));
    lamda(x.dim()-1) = 0.;
    miu(x.dim()-1) = 3.*(y(x.dim()-1) - y(x.dim()-2))/(x(x.dim()-1) - x(x.dim()-2));
    for (int i = 1; i < x.dim()-1; i++) {
        lamda(i) = (x(i) - x(i-1))/(x(i+1) - x(i-1));
        miu(i) = 3.*(1.-lamda(i))*(y(i) - y(i-1))/(x(i) - x(i-1))+3.*lamda(i)*(y(i+1) - y(i))/(x(i+1)-x(i));
    }
    beta(0) = lamda(0)/2.;
    
    for (int i = 1; i < x.dim() - 1; i++) {
        beta(i) = lamda(i)/(2.-beta(i-1)*(1.-lamda(i)));
    }
    alpha(0) = miu(0)/2.;
    for (int i = 1; i < x.dim(); i++) {
        alpha(i) = (miu(i)-(1.-lamda(i))*alpha(i-1))/(2.-(1.-lamda(i))*beta(i-1));
    }
    (*this->item[2])(x.dim()-1) = alpha(x.dim()-1);
    for (int i = x.dim()-2; i >=0; i--) {
        (*this->item[2])(i) = alpha(i) - beta(i)*(*this->item[2])(i+1);
    }
}

template<class T>
const T spline<T>::calculate(const T& x) const{
    int k = -1;
    for (int i = 0; i < (*this->item[0]).dim(); i++) {
        if((x - (*this->item[0])(i)) < 0){
            k = i - 1;
            break;
        }
    }
    T deltax = (*this->item[0])(k+1)-(*this->item[0])(k);
    T alpha1 = (1.+2.*(x - (*this->item[0])(k))/deltax)*(x - (*this->item[0])(k+1))*(x - (*this->item[0])(k+1))/(deltax*deltax);
    T alpha2 = (1.-2.*(x - (*this->item[0])(k+1))/deltax)*(x - (*this->item[0])(k))*(x - (*this->item[0])(k))/(deltax*deltax);
    T beta1 = (x - (*this->item[0])(k))*(x - (*this->item[0])(k+1))*(x - (*this->item[0])(k+1))/(deltax*deltax);
    T beta2 = (x - (*this->item[0])(k+1))*(x - (*this->item[0])(k))*(x - (*this->item[0])(k))/(deltax*deltax);
    if (k==-1){
        return (*this->item[1])((*this->item[1]).dim()-1);
    }
    else {
        return alpha1*(*this->item[1])(k)+alpha2*(*this->item[1])(k+1)+beta1*(*this->item[2])(k)+beta2*(*this->item[2])(k+1);
    }
}


#endif
