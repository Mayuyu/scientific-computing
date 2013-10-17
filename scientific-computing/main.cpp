//
//  main.cpp
//  scientific-computing
//
//  Created by 马 征 on 13-10-17.
//  Copyright (c) 2013年 马 征. All rights reserved.
//

#include <iostream>
#include "gauss_wgts.h"
using namespace std;

int main(int argc, const char * argv[])
{
    dynamicVector<double> a(4,0.),b(4,0.);
    gaulag(a, b, 0.);
    cout<<a<<" "<<endl;
    cout<<b<<endl;
    return 0;
}

