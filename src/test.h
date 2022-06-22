#ifndef _test_h
#define _test_h

#include "Vector/vector_dist.hpp"
#include <iostream>
#include "global.h"
#include "kernel.h"
#include "calculations.h"

void kernel_test(double H);
void output_kernel(double r, double h);
void output_vd(particleset  & vd, int p);
    
#endif