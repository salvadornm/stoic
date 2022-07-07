#ifndef _test_h
#define _test_h

#include "Vector/vector_dist.hpp"
#include <iostream>
#include "global.h"
#include "kernel.h"
#include "calculations.h"

using namespace std;

void kernel_test(double H, Point<3,double> dr);
void output_kernel(double r, double h);
void output_vd(particleset  & vd, int p);
    
#endif