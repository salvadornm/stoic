#ifndef _vector_h
#define _vector_h

#include "Vector/vector_dist.hpp"
#include <iostream>
#include "global.h"
#include <math.h>

double dot_product(vector<double> vector_a, vector<double> vector_b){
    double product = 0;
    for (size_t i = 0; i < 3 ; i++)
    {
        product += vector_a[i] * vector_b[i];
    }
    return product;
}

double dot_product(Point<3,double> vector_a, Point<3,double> vector_b){
    double product = 0;
    for (size_t i = 0; i < 3 ; i++)
    {
        product += vector_a[i] * vector_b[i];
    }
    return product;
}

void cross_product(vector<double> vector_a,vector<double> vector_b, vector<double> & product){
    product = {0,0,0};
    
    product[0] = vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1];
    product[1] = vector_a[2] * vector_b[0] - vector_a[0] * vector_b[2];
    product[2] = vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0];
}

//change coordinates from Normal,Tangential,Tangential
void changeCoordinates(vector<double> xy,vector<double> nt){
    //have x direction vector, have radial vector 
    //θ = cos-1 [ (a · b) / (|a| * |b|) ]

    //double mag_a = sqrt(norm2(xy));
    //double mag_b = sqrt(norm2(nt));
    //double theta = acos((dot_product(xy,nt))/(mag_a * mag_b));
    
    //matrix multiplication (recall T^-1 (i,j) = T(j,i))
}

#endif // _vector_h