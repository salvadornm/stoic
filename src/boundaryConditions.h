#ifndef _boundaryConditions_h
#define _boundaryConditions_h

#include "Vector/vector_dist.hpp"
#include <iostream>
#include "global.h"
#include <math.h>

//determine if point p is inside the cylinder
void topBound(particleset vd, double & vel_z, double & posz, int p, engine eng)
{   
    //(1) Top Wall @ stroke, 
    if (posz - eng.stroke > 0){
    
    //find distance out of bound
    double dz = posz - eng.stroke;
    cout << endl << "top bound posz_1: " << posz;
    posz = eng.stroke - dz;
    cout << " posz_update = " << posz << endl;

    //update velocities
    vel_z = -vel_z;   
    }
}
void pistonBound(particleset vd, double & vel_z, double & posz, int p, engine eng)
{   //is the particle within the z(height) boundaries...
    //assume piston @ BDC and not moving (if moving: vd.getpos(p)[2] < eng.stroke - y) , where y is instantaneous distance from TDC  
    
    //(1) Bottom Wall @ stroke, 
    if (posz < 0){
    //find distance out of bound
    double dz = abs(posz);
    cout << endl << "bottom bound posz_1: " << posz;
    posz = dz;
    cout << " posz_update = " << posz << endl;

    //update velocities
    vel_z = -vel_z;  
    }
}
void sideBC(particleset vd, double & vel_x, double & vel_y, double & posx, double & posy, int p, engine eng, double dt)
{   
    // Get the distance from cylinder wall
    double r_cyl = eng.bore/2;
    Point<2,double> xa {posx,posy};
    Point<2,double> center {r_cyl, r_cyl};
    double r_pos = sqrt(norm2(xa-center));
    double dr_outbound = r_pos - r_cyl;

    // find where the particle hits boundary
    Point<2,double> xb {(posx + vd.getPos(p)[0])/2, (posy + vd.getPos(p)[1])/2}; //position vector of movement
    xb = xb * (r_cyl/sqrt(norm2(xb)));  //scale to radius

    cout << "sideBC- vel x: " << vel_x << " vel y: " << vel_y << endl;
    cout << "sideBC- pos x: " << posx << " pos y: " << posy << endl;

    //if out of bounds in x plane:
    double Ui, Vi, Vr, Ur;
    if (abs(xb[0]) < abs(posx)){
    Ui = vel_y;
    Vi = vel_x;
    Vr = -Vi;
    Ur = Ui - abs(Vi);   // Ui - alpha*Vi

    vel_x = Vr;
    posx += vel_x*dt;
    }

    if (abs(xb[1]) < abs(posy)){
        Vi = vel_y;
        vel_y = -vel_y;
        posy += vel_y*dt;
    }

    cout << "update- vel x: " << vel_x << " vel y: " << vel_y << endl;
    cout << "update- pos x: " << posx << " pos y: " << posy << endl;
/*
    double factor = 2 * ((xb[0]*vel_x + xb[1]*vel_y)/(r_cyl*r_cyl));    
    
    //update velocities
    vel_x = vel_x - factor*xb[0];
    vel_y = vel_y - factor*xb[1];

    posx = posx + vel_x*(dr_outbound/r_pos);
    posy = posy + vel_y*(dr_outbound/r_pos);
    */
}

//return 0 if outside cylinder, 1 if inside cylinder
int inCylinder(particleset vd, double posx, double posy, double posz, int p, engine eng, vector<double> & psi)
{  
    int flag = 1;
    psi = {0.0,0.0,0.0};
    // Get the distance between a and b
    double r_cyl = eng.bore/2;
    Point<2,double> xa {posx,posy};
    Point<2,double> center {r_cyl, r_cyl};

    Point<2,double> dr = xa - center;
    double r2 = norm2(dr);  //norm2 = (sum of the squares)
    double R = sqrt(r2);

    //is the particle within the z(height) boundaries...
    //assume piston @ BDC and not moving (if moving: vd.getpos(p)[2] < eng.stroke - y) , where y is instantaneous distance from TDC  
    if (posz < 0 || posz - eng.stroke > 0){
        flag = 0;   //out of bounds
        cout << "height boundary" << endl;
    }    

    //get point at which particle hits wall
    Point<2,double> xb {(posx + vd.getPos(p)[0])/2, (posy + vd.getPos(p)[1])/2}; //position vector of movement
    xb = xb * (r_cyl/sqrt(norm2(xb)));  //scale to radius
    
    double dz = posz-eng.stroke;
    if (posz < 0){dz = posz;}
    
    psi = {abs(posx)-abs(xb[0]),abs(posy)-abs(xb[1]),dz};

    //is particle within xy plane boundaries....
    if (r_cyl - R < 0) {
        flag = 0;   //out of bounds
        cout << "radial boundary" << endl;
    }
    else{   // if (r_cyl - R > 0)
        flag = 1;   //in bounds
    }

    cout << "point: " << posx << ", " << posy << ", " << posz << endl;
    cout << "wall p: " << xb[0] << ", " << xb[1] << ", " << posz << endl;
    cout << "psi: " << psi[0] << ", " << psi[1] << ", " << psi[2] << endl;
    cout << "r_cyl: " << r_cyl << " R_point: " << R << endl;
    return flag;
}

//determine if point p is inside the cylinder
//return 0 if outside cylinder, 1 if inside cylinder
void initialBoundary(double & posx, double & posy, engine eng, Cfd sim)
{  // Get the distance between a and b
    double r_cyl = eng.bore/2;
    double r = r_cyl * sqrt(((double)rand() / RAND_MAX));
    double theta = ((double)rand()/ RAND_MAX) * 2 * pi;

    posx = r_cyl + r * cos(theta);
    posy = r_cyl + r * sin(theta);
}
#endif // _boundaryConditions_h