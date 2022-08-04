#ifndef _boundaryConditions_h
#define _boundaryConditions_h

#include "Vector/vector_dist.hpp"
#include <iostream>
#include "global.h"
#include <math.h>


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

Point<3,double> wallIntersect(Point <3,double> pos, Point <3,double> pos_new, Point <3,double> center, double r_cyl)
{
    double dx = pos_new[0] - pos[0];
    double dy = pos_new[1] - pos[1];
    double dz = pos_new[2] - pos[2];
    double dxc = pos[0]-center[0];
    double dyc = pos[1]-center[1];
    
    //first find intersection on xy plane
    double a = dx*dx + dy*dy;
    double b = 2*dx*dxc+ 2*dy*dyc;
    double c = dxc*dxc + dyc*dyc - r_cyl*r_cyl;

    double disc = b*b - 4*a*c;
    //double t = (2*c)/(-b + sqrt(disc));
    //double t1 = (-b - sqrt(disc))/(2*a);
    double t2 = (-b + sqrt(disc))/(2*a);

    double x = dx*t2 + pos[0];
    double y = dy*t2 + pos[1];

    //then find at what height the line is at that xy position
    double z = dz*t2 + pos[2];
    Point<3,double> pos_wall{x,y,z};
    return pos_wall;
}

//return 0 if outside cylinder, 1 if inside cylinder
int inCylinder(particleset vd, Point <3,double> pos, Point <3,double> pos_new, int p, engine eng, Point<3,double> & psi)
{  
    int flag = 1;
    psi = {0.0,0.0,0.0};

    // set up geometries
    double r_cyl = eng.bore/2;
    Point<3,double> cyl_center {r_cyl, r_cyl, pos_new[2]};
    double r_pos = sqrt(norm2(pos_new - cyl_center));

    Point<3,double> pos_wall {0,0,0};

    //is the particle within the z(height) boundaries...
    //assume piston @ BDC and not moving (if moving: vd.getpos(p)[2] < eng.stroke - y) , where y is instantaneous distance from TDC  
    if (pos_new[2] < 0 || pos_new[2] - eng.stroke > 0){
        flag = 0;   //out of bounds
        cout << "height boundary" << endl;
    }        
    
    //is particle within xy plane boundaries....
    if (r_cyl - r_pos < 0) {
        flag = 0;   //out of bounds
        cout << "radial boundary" << endl;
        //get point at which particle hits wall
        pos_wall = wallIntersect(pos, pos_new, cyl_center, r_cyl);
    }

    psi = pos_new - pos_wall;

    //output for testing
    cout << "point: " << pos_new[0] << ", " << pos_new[1] << ", " << pos_new[2] << endl;
    cout << "wall p: " << pos_wall[0] << ", " << pos_wall[1] << ", " << pos_wall[2] << endl;
    cout << "psi: " << psi[0] << ", " << psi[1] << ", " << psi[2] << endl;
    cout << "r_cyl: " << r_cyl << " R_point: " << r_pos << endl;
    return flag;
}

//determine if point p is inside the cylinder <-- looks like this is working!
void topBound(Point <3,double> & pos_new, vector <double> &vel, engine eng)
{   
    //(1) Top Wall @ stroke, 
    if (pos_new[2] - eng.stroke > 0){
    
    //find distance out of bound
    double dz = pos_new[2] - eng.stroke;
    cout << endl << "top bound posz_1: " << pos_new[2];
    pos_new[2] = eng.stroke - dz;
    cout << " posz_update = " << pos_new[2] << endl;

    //update velocities
    vel[2] = -vel[2];   
    }
}


void pistonBound(Point <3,double> & pos_new, vector <double> &vel, engine eng)
{   //assume piston @ BDC and not moving (if moving: vd.getpos(p)[2] < eng.stroke - y) , where y is instantaneous distance from TDC  
    
    //(1) Bottom Wall @ stroke, 
    if (pos_new[2] < 0){
    //find distance out of bound
    double dz = abs(pos_new[2]);
    cout << endl << "bottom bound posz_1: " << pos_new[2];
    pos_new[2] = dz;
    cout << " posz_update = " << pos_new[2] << endl;

    //update velocities
    vel[2] = -vel[2];  
    }
}


void sideBC(particleset vd, Point <3,double> & pos_zero, Point <3,double> & pos_new, vector <double> &vel, int p, engine eng, double dt)
{   
    // set up geometries
    double r_cyl = eng.bore/2;
    Point<3,double> cyl_center {r_cyl, r_cyl, pos_new[2]};
    double r_pos = sqrt(norm2(pos_new - cyl_center));

    // Get the distance from cylinder wall
    Point<3,double> pos_wall = wallIntersect(pos_zero, pos_new, cyl_center, r_cyl);

    cout << "sideBC- vel x: " << vel[0] << " vel y: " << vel[1] << " vel z: " << vel[2] << endl;
    cout << "sideBC- pos x: " << pos_wall[0] << " pos y: " << pos_wall[1] << " pos z: " << pos_wall[2] << endl;

    // projection of v onto n
    double factor = 2 * ((pos_wall[0]*vel[0] + pos_wall[1]*vel[1])/(r_cyl*r_cyl));    
    
    //update velocities
    vel[0] = vel[0] - factor*pos_wall[0];
    vel[1] = vel[1] - factor*pos_wall[1];

    pos_new[0] = pos_new[0] + vel[0]*dt;    //DONT THINK THIS IS RIGHT
    pos_new[1] = pos_new[2] + vel[1]*dt;
    //bounce doesnt impact height

/*
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
*/

    cout << "update- vel x: " << vel[0] << " vel y: " << vel[1] << endl;
    cout << "update- pos x: " << pos_wall[0] << " pos y: " << pos_wall[1] << endl;
}

#endif // _boundaryConditions_h