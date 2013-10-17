//
//  complex.h
//  scientific-computing
//
//  Created by 马 征 on 13-10-17.
//  Copyright (c) 2013年 马 征. All rights reserved.
//

#ifndef scientific_computing_complex_h
#define scientific_computing_complex_h

#include <iostream>

using namespace std;

class complex {
    double real;
    double image;
public:
    complex(double x = 0, double y = 0):real(x),image(y){}
    complex(const complex & T):real(T.real),image(T.image){}
    const double re()const{return real;}
    const double im()const{return image;}
    const complex& operator=(const complex&);
    const complex& operator+=(const complex&);
    const complex& operator-=(const complex&);
    const complex& operator*=(const complex&);
    const complex& operator/=(const double&);
    const complex& operator/=(const complex&);
    ~complex(){}
};

const complex& complex::operator=(const complex& c){
    real = c.real;
    image = c.image;
    return *this;
}

const complex& complex::operator+=(const complex& c){
    real += c.real;
    image += c.image;
    return *this;
}

const complex& complex::operator-=(const complex& c){
    real -=c.real;
    image -=c.image;
    return *this;
}

const complex& complex::operator*=(const complex& c){
    double temp = real;
    real = real*c.real - image*c.image;
    image = temp*c.image + image*c.real;
    return *this;
}

const complex& complex::operator/=(const double& d){
    real /= d;
    image /= d;
    return *this;
}

double abs2(const complex& c){
    return c.re()*c.re()+c.im()*c.im();
}

complex operator!(const complex& c){
    return complex(c.re(), -c.im());
}

const complex& complex::operator/=(const complex& c){
    return *this *= (!c) /= abs2(c);
}

const complex operator+(const complex& a, const complex& b){
    return complex(a)+=b;
}

const complex operator-(const complex& c){
    return complex(-c.re(), -c.im());
}

const complex operator-(const complex& a, const complex& b){
    return complex(a)-=b;
}

const complex operator*(const complex& a, const complex& b){
    return complex (a) *= b;
}

const complex operator/(const complex&a, const complex& b){
    return complex(a) /= b;
}

ostream& operator <<(ostream& out, const complex& c){
    out << c.re() << "+" << c.im() << "i";
    return out;
}


#endif
