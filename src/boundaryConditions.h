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
    //double t1 = (-b - sqrt(disc))/(2*a);
    double t2 = (-b + sqrt(disc))/(2*a);

    double x = dx*t2 + pos[0];
    double y = dy*t2 + pos[1];

    //then find at what height the line is at that xy position
    double z = dz*t2 + pos[2];
    Point<3,double> pos_wall{x,y,z};
    return pos_wall;
}

void hitFirst(vector<double> vel, Point <3,double> pos_zero, Point <3,double> pos_wall, Point<3,double> & psi)
{
    double deltaT;
    for (size_t i = 0; i < 3 ; i++) 
    {  
        deltaT = (pos_wall[i] - pos_zero[i])/vel[i];
        psi[i] = deltaT;
    }
}

//return 0 if outside cylinder, 1 if inside cylinder
int inCylinder(vector<double> vel, Point <3,double> pos, Point <3,double> pos_new, engine eng, Point<3,double> & psi)
{  
    int flag = 1;
    psi = {0.0,0.0,0.0};

    // set up geometries
    double r_cyl = eng.bore/2;
    Point<3,double> cyl_center {r_cyl, r_cyl, pos_new[2]};
    double r_pos = sqrt(norm2(pos_new - cyl_center));

    Point<3,double> pos_wall {0,0,0};
    pos_wall = wallIntersect(pos, pos_new, cyl_center, r_cyl);

    //is particle within xy plane boundaries....
    if (r_cyl - r_pos < 0) {
        flag = 0;   //out of bounds
        //cout << "radial boundary" << endl;
    }

    //is the particle within the z(height) boundaries...
    //if moving: vd.getpos(p)[2] < eng.stroke - y) , where y is instantaneous distance from TDC  
    if (pos_new[2] < eng.s_inst){
        flag = 0;   //out of bounds
    }    
    else if(pos_new[2] - eng.height > 0){
        flag = 0;   //out of bounds
    }

    if (vel[2] > 0){pos_wall[2] = eng.height; }
    else {pos_wall[2] = eng.s_inst; }

    hitFirst(vel, pos, pos_wall, psi);

    //output for testing
    //output_bc_props(vel, pos, pos_new, pos_wall, eng, psi);
    //cout << "r_cyl: " << r_cyl << " R_point: " << r_pos << endl;
    
    return flag;
}

//return 0 if outside cylinder, 1 if inside cylinder USED SPECIFICALLY FOR PISTON MOVEMENT
int inCylinder(Point <3,double> pos, engine eng)
{  
    int flag = 1;
    //is the particle within the z(height) boundaries... 
    if (pos[2] < eng.s_inst){
        flag = 0;   //out of bounds
    }    
    else if(pos[2] - eng.height > 0){
        flag = 0;   //out of bounds
    }

    return flag;
}

//determine if point p is inside the cylinder
void topBound(Point <3,double> & pos_new, vector <double> &vel, engine eng)
{   
    //(1) Top Wall @ height, 
    if (pos_new[2] - eng.height > 0){
    
    //find distance out of bound
    double dz = pos_new[2] - eng.height;
    pos_new[2] = eng.height - dz;

    //update velocities
    vel[2] = -vel[2];   
    }
}

void pistonBound(Point <3,double> & pos_new, vector <double> &vel, engine eng)
{   // assume piston @ BDC, and piston geometry is flat
    // if moving: vd.getpos(p)[2] < eng.stroke - y, where y is instantaneous distance from TDC  
    
    //(1) Bottom Wall @ stroke, 
    if (pos_new[2] < eng.s_inst){
    //find distance out of bound
    double dz = abs(pos_new[2] - eng.s_inst);
    pos_new[2] = eng.s_inst + dz;

    //update normal velocity to match piston speed
    vel[2] = eng.smp;  
    }
}


void sideBC(particleset vd, Point <3,double> & pos_zero, Point <3,double> & pos_new, vector <double> &vel, int p, engine eng, double dt)
{   
    // set up geometries
    double r_cyl = eng.bore/2;
    Point<3,double> cyl_center {r_cyl, r_cyl, pos_new[2]};
    double r_pos = sqrt(norm2(pos_new - cyl_center));
    
    if (r_pos > r_cyl){
        
    // Get the distance from cylinder wall
    Point<3,double> pos_wall = wallIntersect(pos_zero, pos_new, cyl_center, r_cyl);

    //make changing from x,y,z to norm, tan1, tan2 coordinate systems a function
    // projection of v onto n
    Point<3,double> newV = pos_new - pos_wall;
    Point<3,double> radialV = pos_wall - cyl_center;

    double factor = 2 * ((radialV[0]*vel[0] + radialV[1]*vel[1])/(r_cyl*r_cyl));

    //update velocities
    vel[0] = vel[0] - factor*radialV[0];
    vel[1] = vel[1] - factor*radialV[1];

    //update position (using different velocities)
    double factorP = 2 * ((radialV[0]*newV[0] + radialV[1]*newV[1])/(r_cyl*r_cyl));  //dot product  
    newV[0] = newV[0] - factorP*radialV[0];
    newV[1] = newV[1] - factorP*radialV[1];
    pos_new[0] = pos_wall[0] + newV[0]*dt;   
    pos_new[1] = pos_wall[1] + newV[1]*dt;

    //bounce doesnt impact height

    //cout << "update- pos x: " << pos_new[0] << " pos y: " << pos_new[1] << endl;
    //cout << "update- vel x: " << vel[0] << " vel y: " << vel[1] << " vel z: " << vel[2] << endl;
    }
}

#endif // _boundaryConditions_h